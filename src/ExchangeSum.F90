!Routine for constructing the a/c lattice sum list, and the orbital pairs that need to go with the
!exchange interaction.
!IN: Cell info
!    dCoDensTol - Tolerance for inclusion of codensity in list, based on |(ij|ij)|^0.5
!OUT: Ex_Codens, which contains:
!   - Number of unique translation vectors between all groups in codensity: nTransVecs
!   - List of these translation vectors (vector of 3 integers)  TransVecs(:,:)
!   - Number of pairs of groups which contribute
!   - Number of group pairs where first group index is in the unit cell (these are indexed first)
!   - Vector of orbital pairs, containing:
!       - Indices of groups (supercell labeled)
!       - Supercell translation vector between the two groups (c)
!       - Max value of |(1 2_c | 1 2_c)|^0.5 between functions in each group   (used for CS screening)
!       - Approximate max radii of codensity functions between groups (used for b lattice sum screening)
!       - Approximate center of codensity between groups (used for b lattice sum screening)

!Two parameters - whether to choose the codensity list by the max range of the contributing shells (tMakeListByRange),
!...and if so, whether to nevertheless calculate and store their exchange integrals (tCreateCSList).

!Locally, also output largest translation vector in terms of sum of squared distance
subroutine ConstructLattSumOrbPairs(Ex_Codens,Lattice,UnitCell,Supercell,dCoDensTol)
    use CvsFDatatypes
    use error_mod, only: stop_all
    use FortCons, only: dp
    implicit none
    !Arguments
    type(exchange_densities_t), intent(out) :: Ex_Codens
    type(lattice_t) , intent(in) :: Lattice
    type(unit_cell_t), intent(in)  :: UnitCell
    type(supercell_t), intent(in) :: Supercell
    real(dp), intent(in) :: dCoDensTol
    !Unpacked arrays
    type(GroupInfo), pointer :: Supercell_groups(:)
    type(GroupInfo), pointer :: Unitcell_groups(:)
    real(dp), pointer :: Supercell_data(:)
    real(dp), pointer :: Supercell_centers(:,:),Unitcell_coords(:,:)
    !local arguments
    integer :: iMaxCodensTrans  !maximum translation in any one direction from gamma-point to consider
    integer :: UnitCellFns,SupercellFns,abf,bbf,nAtoms,c
    integer :: ngroups_super,ngroups_unit,nTrans,iNumGroups
    integer :: ic,nlargestshell,strides(4),iMaxTrans,cbf,dbf
    integer :: i,a,b,nbf_a,nbf_b,lx,lxd,x,ly,lyd,y,z,lz,lzd
    real(dp) :: L,L_x,L_y,L_z,Rc,zTrans(3),TransVec(3),magint,maxrange,mindim
    real(dp) :: rangefuncs,distance,center_a(3),center_b(3),maxint
    integer, allocatable :: MaxTransVecs(:,:)
    real(dp), allocatable :: int_temp(:,:,:,:)
    logical ::  tIncTrans
    !parameters
    character(len=*), parameter :: t_r='ConstructLattSumOrbPairs'
    logical, parameter :: tMakeListByRange = .true.  
    logical, parameter :: tCreateCSList = .true.
    !external functions
    integer :: create_integral_context
    
    UnitCellFns = UnitCell%OrbBasis%nFn
    SupercellFns = Supercell%OrbBasis%nFn
    
    !Extract the basis set over the supercell
    call c_f_pointer(Supercell%OrbBasis%Groups%pData, Supercell_groups, [Supercell%OrbBasis%Groups%nSize])
    !Extract the basis set over the unit cell 
    call c_f_pointer(Unitcell%OrbBasis%Groups%pData, Unitcell_groups, [Unitcell%OrbBasis%Groups%nSize])
    nAtoms = Unitcell%Coords%nSize
    call c_f_pointer(Unitcell%Coords%pData, Unitcell_coords, [3,Unitcell%Coords%nSize])
    !Extract the data for the supercell functions (You never get this directly, but rather it is pointed to from group objects)
    !When pointing to elements in this array, beware that it points to the element before it, therefore, you need to add 1!1
    call c_f_pointer(Supercell%OrbBasis%Data%pData, Supercell_data, [Supercell%OrbBasis%Data%nSize])
    call c_f_pointer(Supercell%OrbBasis%Centers%pData, Supercell_centers, [3,Supercell%OrbBasis%Centers%nSize])

    ngroups_super = Supercell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in supercell: ",ngroups_super
    write(6,*) "Number of functions in supercell: ",SupercellFns
    ngroups_unit = Unitcell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in unitcell: ",ngroups_unit

    write(6,"(A,3F12.6)") "SC Lattice_x = ",Supercell%T(1,1),Supercell%T(2,1),Supercell%T(3,1)
    write(6,"(A,3F12.6)") "SC Lattice_y = ",Supercell%T(1,2),Supercell%T(2,2),Supercell%T(3,2)
    write(6,"(A,3F12.6)") "SC Lattice_z = ",Supercell%T(1,3),Supercell%T(2,3),Supercell%T(3,3)
    write(6,*) ""
    write(6,*) "Atom coordinates: "
    do i=1,nAtoms
        write(6,"(A,I5,3F12.6)") "Atom ",i,Unitcell_coords(:,i)
    enddo
    
!In future, there should be a routine here for constructing the optimal Rc 
    !Want Rc .le. L/2
    !Calculate average L for supercell
    L_x = sqrt(Supercell%T(1,1)**2 + Supercell%T(2,1)**2 + Supercell%T(3,1)**2)
    L_y = sqrt(Supercell%T(1,2)**2 + Supercell%T(2,2)**2 + Supercell%T(3,2)**2)
    L_z = sqrt(Supercell%T(1,3)**2 + Supercell%T(2,3)**2 + Supercell%T(3,3)**2)
    L = (L_x + L_y + L_z) / 3.0_dp
    write(6,*) "Average length of supercell dimensions is: ",L
    Rc = L/2.0_dp 
    write(6,*) "Truncated coulomb potential kernel cutoff: ",Rc

    ic = create_integral_context(0,0,1.0e-10_dp)
    !Plain coulomb
    !call assign_integral_kernel(ic,3,0,0)  !3 is for coulomb kernal
    !Truncated coulomb
    call assign_integral_kernel(ic,6,0,Rc)  !6 is for truncated coulomb kernal
    write(6,"(A)") "Constructing screened list of codensities for exchange summation"
    if(tMakeListbyRange) then
        write(6,"(A)") "This screening be based on a range criteria of the orbitals which make up the codensity"
        if(tCreateCSList) write(6,"(A)") "The Cauchy-Schwartz criteria will also be saved for use later on"
    else
        write(6,"(A)") "This screening will be based on a Cauchy-Schwarz criterion"
        write(6,"(A)") "This will be slower than the range based criteria"
    endif
    
    !Find an approximate guess for the number of supercell lattice vectors to use in each direction
    !A highly conservative guess will be int(largest function range*2/smallest supercell dimension) + 1
    !This should be a good approximation to the greatest distance we need to go out
    maxrange = 0.0_dp
    do i = 1,ngroups_super
        !write(6,*) Supercell_groups(i)%iRange,Supercell_data(Supercell_groups(i)%iRange+1)
        if(Supercell_data(Supercell_groups(i)%iRange+1).gt.maxrange) then
            maxrange = Supercell_data(Supercell_groups(i)%iRange+1)
        !    write(6,*) "MaxRange: ",Supercell_data(Supercell_groups(i)%iRange+1)
        endif
    enddo
    if(maxrange.lt.1.0e-5_dp) call stop_all(t_r,'Error here')
    write(6,*) "Largest range of any one basis function in the supercell is approximately: ",maxrange
    mindim = L_x 
    if(L_y.lt.mindim) then
        mindim = L_y
    endif
    if(L_z.lt.mindim) then
        mindim = L_z
    endif
    write(6,*) "Smallest supercell dimension is: ",mindim
    iMaxCodensTrans = int(maxrange*2.0_dp/mindim) + 1
    write(6,*) "Considering all supercell translation vectors in each direction up to ",iMaxCodensTrans

    if((.not.tMakeListbyRange).or.tCreateCSList) then
        !What is the largest number of integrals (ij|ij) that can come out of two shells?
        !This is so we can store them when we calculate the exchange integrals.
        nlargestshell = 0
        do i = 1,ngroups_super
            if(Supercell_groups(i)%nFn.gt.nlargestshell) nlargestshell = Supercell_groups(i)%nFn
        enddo
        write(6,*) "Largest number of functions in a group: ",nlargestshell
        !How many integrals does this mean that we can get out of an (ij|ij) set?
        allocate(int_temp(nlargestshell,nlargestshell,nlargestshell,nlargestshell)) !This is to temporarily hold the integrals
        strides = (/ 1 , nlargestshell , nlargestshell**2 , nlargestshell**3 /)
    endif

!Allocate memory to hold maximum possible translation vectors that we need to consider
!iMaxCodensTrans gives the number of translations in a given direction. Find the maximum number of translation vectors to consider
    iMaxTrans = (2*iMaxCodensTrans+1)**3
    write(6,*) "Maximum number of translation vectors to search through in constructing group list: ",iMaxTrans
    allocate(MaxTransVecs(3,iMaxTrans))
    MaxTransVecs(:,:) = 0
    zTrans(:) = 0.0_dp  !The null translation vector
    nTrans = 0  !Index of current translation to consider
    iNumGroups = 0  !Number of pairs of groups we are including

    !First, count the number of group pairs we want, and construct the temporary list of translation vectors
    !loop over all translation vectors
    do lxd = 1,2    !Which direction to go in
        do lx = 0,iMaxCodensTrans
            if((lx.eq.0).and.(lxd.eq.2)) cycle
            if(lxd.eq.1) then
                x = lx
            else
                x = -lx
            endif
            do lyd = 1,2
                do ly = 0,iMaxCodensTrans
                    if(lyd.eq.1) then
                        y = ly
                    else
                        y = -ly
                    endif
                    if((ly.eq.0).and.(lyd.eq.2)) cycle
                    do lzd = 1,2
                        do lz = 0,iMaxCodensTrans
                            if((lz.eq.0).and.(lzd.eq.2)) cycle
                            if(lzd.eq.1) then
                                z = lz
                            else
                                z = -lz
                            endif

                            if((.not.tIncTrans).and.(z.ne.0).and.(z.ne.-1)) then
                                !We didn't want any of the group pairs in the last translation,
                                !And we have only gone out further now. Cycle through the rest of this
                                !lz set...
                                !This shouldn't result in any less pairs added to the list, and just be
                                !and efficiency improvement
                                cycle
                            endif

                            !Calculate translation vector for groups b
                            TransVec(:) = 0.0_dp    !3-vector of displacements for group b
                            TransVec(:) = x*Supercell%T(:,1) + y*Supercell%T(:,2) + z*Supercell%T(:,3)
                            tIncTrans = .false. !flag to indicate whether to include this translation in the list
                            !write(6,"(A,3I5)") "Considering SC translation: ",x,y,z

                            do a = 1,ngroups_super
                                !write(6,*) "Considering gamma-point group: ",a
                                nbf_a = Supercell_groups(a)%nFn

                                !Now loop over functions in this translated supercell
                                do b = 1,ngroups_super
                                    nbf_b = Supercell_groups(b)%nFn

                                    if(tMakeListbyRange) then
                                        !if |center_a - (center_b + T_c)| < Range_a + Range_b, then accept
                                        center_a(1:3) = Supercell_centers(:,Supercell_groups(a)%iCen+1)
                                        center_b(1:3) = Supercell_centers(:,Supercell_groups(b)%iCen+1)+TransVec(:)
                                        distance = (center_a(1)-center_b(1))**2
                                        distance = distance + (center_a(2)-center_b(2))**2
                                        distance = distance + (center_a(3)-center_b(3))**2
                                        distance = sqrt(distance)
                                        !write(6,*) "Center a: ",center_a(:)
                                        !write(6,*) "Center b: ",center_b(:)
                                        !write(6,*) "Distance between functions: ",a,b,distance
                                        rangefuncs = Supercell_data(Supercell_groups(a)%iRange+1) + Supercell_data(Supercell_groups(b)%iRange+1)

                                        if(distance.lt.rangefuncs) then
                                            !include these groups in list
                                            if(.not.tIncTrans) then
                                                !We want to include this translation vector
                                                nTrans = nTrans + 1 !Increment list of lattice vectors
                                                MaxTransVecs(1,nTrans) = x
                                                MaxTransVecs(2,nTrans) = y
                                                MaxTransVecs(3,nTrans) = z
                                                tIncTrans = .true.
                                            endif
                                            iNumGroups = iNumGroups + 1
                                        endif
                                    else
                                        !Include based on the max value of |(ab^c|ab^c)|^0.5
                                        int_temp(:,:,:,:) = 0.0_dp
                                        !Calculate set of integrals between all functions between these groups
                                        call eval_group_int2e_tra_incr(int_temp(1,1,1,1),strides(1),    &
                                            1.0_dp,a-1,zTrans(1),b-1,TransVec(1),a-1,zTrans(1),b-1,     &
                                            TransVec(1),Supercell%OrbBasis,ic)

                                        !Now, we want to run through the diagonal exchange integrals, to check for
                                        !whether to include them or not
                                        loop: do abf = 1,nbf_a
                                            do bbf = 1,nbf_b
                                                magint = sqrt(abs(int_temp(abf,bbf,abf,bbf)))
                                                if(magint.gt.dCoDensTol) then
                                                    !We want to include consideration of this group
                                                    if(.not.tIncTrans) then
                                                        !We want to include this translation vector
                                                        nTrans = nTrans + 1 !Increment list of lattice vectors
                                                        MaxTransVecs(1,nTrans) = x
                                                        MaxTransVecs(2,nTrans) = y
                                                        MaxTransVecs(3,nTrans) = z
                                                        tIncTrans = .true.
                                                    endif
                                                    iNumGroups = iNumGroups + 1
                                                    !write(16,*) x**2+y**2+z**2,magint
                                                    exit loop
                                                endif
                                            enddo
                                        enddo loop
                                    endif

                                enddo   !end b
                            enddo   ! end a
                        enddo   !lz
                    enddo   !lzd
                enddo   !ly
            enddo   !lyd
        enddo   !lx
    enddo   !lxd

    write(6,*) "SC Translation vectors that we need to consider for codensities: "
    do i = 1,nTrans
        write(6,"(3I6)") MaxTransVecs(:,i)
    enddo
    write(6,"(A,I9,A,I9,A)") "This corresponds to ",nTrans," SC translations out of a possible searched ",iMaxTrans," translations."
    write(6,*) ""
    write(6,"(A,I9,A,I9)") "Total number of group pairs to consider for each electron in exchange sum: ",iNumGroups,    &
        " out of a possible ",(ngroups_super**2)*iMaxTrans

    !Now, move this list into the appropriate derived type
    Ex_Codens%nTransVecs = nTrans
    Ex_Codens%ngroup_pairs = iNumGroups
    Ex_Codens%ngroup_pairs_UC = 0 
    allocate(Ex_Codens%iTransVecs(3,nTrans))
    Ex_Codens%iTransVecs(:,:) = 0
    do i = 1,nTrans
        Ex_Codens%iTransVecs(:,i) = MaxTransVecs(:,i)
    enddo
    deallocate(MaxTransVecs)

    !Now allocate the group pair list
    allocate(Ex_Codens%grouppairlist%group_labels(2,iNumGroups))
    Ex_Codens%grouppairlist%group_labels(:,:) = 0
    allocate(Ex_Codens%grouppairlist%transVec(3,iNumGroups))
    Ex_Codens%grouppairlist%transVec(:,:) = 0.0_dp
    if((tMakeListbyRange.and.tCreateCSList).or.(.not.tMakeListbyRange)) then
        !We are not using the CS criteria to select the groups, but store them anyway
        allocate(Ex_Codens%grouppairlist%CS_Fac(iNumGroups))
        Ex_Codens%grouppairlist%CS_Fac(:) = 0.0_dp
    endif
    allocate(Ex_Codens%grouppairlist%den_radius(iNumGroups))
    Ex_Codens%grouppairlist%den_radius(:) = 0.0_dp
    allocate(Ex_Codens%grouppairlist%den_center(3,iNumGroups))
    Ex_Codens%grouppairlist%den_center(:,:) = 0.0_dp

    iNumGroups = 0  !Reset this counter
    !Now fill these arrays!
    do a = 1,ngroups_super
        !write(6,*) "Considering gamma-point group: ",a
        nbf_a = Supercell_groups(a)%nFn
        if(a.eq.ngroups_unit+1) then
            !We have just started moving into functions outside the unit cell
            Ex_Codens%ngroup_pairs_UC = iNumGroups 
        endif

        !Now loop over functions in this translated supercell
        do b = 1,ngroups_super
            nbf_b = Supercell_groups(b)%nFn

            do c = 1,nTrans !Loop over translations of b
                TransVec(:) = Ex_Codens%iTransVecs(1,c)*Supercell%T(:,1) +   &
                              Ex_Codens%iTransVecs(2,c)*Supercell%T(:,2) +   &
                              Ex_Codens%iTransVecs(3,c)*Supercell%T(:,3)

                !Now determine if we want to store this group pair
                if(tMakeListbyRange) then
                    !if |center_a - (center_b + T_c)| < Range_a + Range_b, then accept
                    center_a(1:3) = Supercell_centers(:,Supercell_groups(a)%iCen+1)
                    center_b(1:3) = Supercell_centers(:,Supercell_groups(b)%iCen+1)+TransVec(:)
                    distance = (center_a(1)-center_b(1))**2
                    distance = distance + (center_a(2)-center_b(2))**2
                    distance = distance + (center_a(3)-center_b(3))**2
                    distance = sqrt(distance)
                    !write(6,*) "Center a: ",center_a(:)
                    !write(6,*) "Center b: ",center_b(:)
                    !write(6,*) "Distance between functions: ",a,b,distance
                    rangefuncs = Supercell_data(Supercell_groups(a)%iRange+1) + Supercell_data(Supercell_groups(b)%iRange+1)

                    if(distance.lt.rangefuncs) then
                        !include these groups in list
                        iNumGroups = iNumGroups + 1
                        if(iNumGroups.gt.Ex_Codens%ngroup_pairs) call stop_all(t_r,'Error counting')

                        Ex_Codens%grouppairlist%group_labels(1,iNumGroups) = a
                        Ex_Codens%grouppairlist%group_labels(2,iNumGroups) = b
                        Ex_Codens%grouppairlist%transVec(:,iNumGroups) = TransVec(:)
                        if(tCreateCSList) then
                            !Find max CS integral between the groups
                            int_temp(:,:,:,:) = 0.0_dp
                            !Calculate set of integrals between all functions between these groups
                            call eval_group_int2e_tra_incr(int_temp(1,1,1,1),strides(1),    &
                                1.0_dp,a-1,zTrans(1),b-1,TransVec(1),a-1,zTrans(1),b-1,     &
                                TransVec(1),Supercell%OrbBasis,ic)

                            !Now, we want to run through the diagonal exchange integrals, to check for
                            !whether to include them or not
                            maxint = 0.0_dp
                            do abf = 1,nbf_a
                                do bbf = 1,nbf_b
                                    !We have to just go through all the integrals here, since 
                                    !we want to store the maximum one
                                    magint = sqrt(abs(int_temp(abf,bbf,abf,bbf)))
                                    if(magint.gt.maxint) maxint = magint
                                enddo
                            enddo
                            Ex_Codens%grouppairlist%CS_Fac(iNumGroups) = maxint
                        endif
                        !Find approximate radius of this codensity TODO
                        !Simplest to start with - midpoint of functions, and radius = distance between the functions
                    endif
                else
                    !Include based on the max value of |(ab^c|ab^c)|^0.5
                    int_temp(:,:,:,:) = 0.0_dp
                    !Calculate set of integrals between all functions between these groups
                    call eval_group_int2e_tra_incr(int_temp(1,1,1,1),strides(1),    &
                        1.0_dp,a-1,zTrans(1),b-1,TransVec(1),a-1,zTrans(1),b-1,     &
                        TransVec(1),Supercell%OrbBasis,ic)

                    !Now, we want to run through the diagonal exchange integrals, to check for
                    !whether to include them or not
                    maxint = 0.0_dp
                    do abf = 1,nbf_a
                        do bbf = 1,nbf_b
                            !We have to just go through all the integrals here, since 
                            !we want to store the maximum one
                            magint = sqrt(abs(int_temp(abf,bbf,abf,bbf)))
                            if(magint.gt.maxint) maxint = magint
                        enddo
                    enddo
                    if(maxint.gt.dCoDensTol) then
                        !We want to include consideration of this group
                        iNumGroups = iNumGroups + 1
                        if(iNumGroups.gt.Ex_Codens%ngroup_pairs) call stop_all(t_r,'Error counting')
                        !write(16,*) x**2+y**2+z**2,magint
                        Ex_Codens%grouppairlist%group_labels(1,iNumGroups) = a
                        Ex_Codens%grouppairlist%group_labels(2,iNumGroups) = b
                        Ex_Codens%grouppairlist%transVec(:,iNumGroups) = TransVec(:)
                        Ex_Codens%grouppairlist%CS_Fac(iNumGroups) = maxint
                        !Now also calculate radius and center of codensity!
                    endif

                endif

            enddo   !Lattice sum
        enddo   !b
    enddo   !a

    if(iNumGroups.ne.Ex_Codens%ngroup_pairs) call stop_all(t_r,'Did not store expected number of pairs')

    write(6,*) "Finished storing codensity list for exchange calculation..."

end subroutine ConstructLattSumOrbPairs

subroutine ExchangeSum(ExEnergy,Exchange,Lattice,UnitCell,Supercell,Density)
    use CvsFDatatypes
    use error_mod, only: stop_all
    use FortCons, only: dp
    implicit none
    !Arguments
    type(op_matrix_t), intent(out) :: Exchange
    real(dp), intent(out) :: ExEnergy
    type(op_matrix_t), intent(in) :: Density
    type(lattice_t) , intent(in) :: Lattice
    type(unit_cell_t), intent(in)  :: UnitCell
    type(supercell_t), intent(in) :: Supercell

    type(Exchange_densities_t) :: Ex_Codens

    !Unpacked arguments
    real(dp), pointer :: DenMat(:,:)
    real(dp), pointer :: ExMat(:,:)
    type(GroupInfo), pointer :: Supercell_groups(:)
    type(GroupInfo), pointer :: Unitcell_groups(:)

    !Local variables
    real(dp) :: dDenDim   !Dimension of density matrix
    real(dp) :: mu_nu_CS
    real(dp) :: Rc,L_x,L_y,L_z,L,TransVec(3),BplusCTrans(3),mu_nu_trans(3),lam_sig_trans(3),zTrans(3)
    integer :: UnitCellFns,lam_ind,sig_ind,nu_ind,mu_ind,lam_sig,lambf,mubf,nubf,sigbf
    integer :: lx,ly,lz,lxd,lyd,lzd,mu_nu,nbf_lam,nbf_mu,nbf_nu,nbf_sig
    integer :: SupercellFns,ngroups_super,ngroups_unit,nlargestshell,x,y,z
    integer :: lam,mu,nu,sig,i,ic,ierr,iMaxCodensTrans 
    real(dp), allocatable :: UnpackedDM(:,:)    !SS x SS unpacked DM
    real(dp), allocatable :: Int_Temp(:,:,:,:)
    integer, allocatable :: FunctionInd_Shell(:)
    integer :: strides(4)
    logical :: tIncTrans,tLattSum
    character(len=*), parameter :: t_r='ExchangeSum'
    real(dp) , parameter :: CS_CodensThresh = 1.0e-7_dp
    real(dp) , parameter :: NF_ScreenTol = 1.0e-8_dp
    real(dp) , parameter :: DM_ScreenTol = 1.0e-8_dp
    logical , parameter :: Far_fieldScreen = .true.

    !external functions
    integer :: create_integral_context
    real(dp) :: DDOT

    UnitCellFns = UnitCell%OrbBasis%nFn
    !write(6,*) "Number of basis functions in unit cell: ",UnitCellFns
    SupercellFns = Supercell%OrbBasis%nFn
    !write(6,*) "Number of basis functions in supercell: ",SupercellFns

    !write(6,*) "Number of unit cells in supercell: ",Supercell%Size(:)
    if((SupercellFns-(Supercell%Size(1)*Supercell%Size(2)*  &
        Supercell%Size(3))*UnitCellFns).ne.0) then
        write(6,*) "SupercellFns: ",SupercellFns
        write(6,*) "prod(Supercell%Size(:))*UnitCellFns: ",Supercell%Size(1)*   &
            Supercell%Size(2)*Supercell%Size(3)*UnitCellFns
        call stop_all(t_r,"Error in expected number of basis functions in supercell")
    endif

    !Extract size of passed in density
    if(Density%nrows.ne.UnitCellFns) then
        call stop_all(t_r,'First dim of density matrix not of expected size')
    endif
    if(Density%ncols.ne.SupercellFns) then
        call stop_all(t_r,'Second dim of density matrix not of expected size')
    endif
    if(Density%Data%nSize.ne.(UnitCellFns*SupercellFns)) then
        call stop_all(t_r,"Unexpected size of density matrix")
    endif
    if(Density%Data%nSize.ne.Exchange%Data%nSize) then
        call stop_all(t_r,"Density matrix not of same size as exchange matrix")
    endif
    !Now extract density matrix element
    call c_f_pointer(Density%Data%pData, DenMat, [UnitCellFns,SupercellFns])
    call c_f_pointer(Exchange%Data%pData, ExMat, [UnitCellFns,SupercellFns])

    !write(6,*) "Size of matrices as expected"

    !Unpack symmetry
    allocate(UnpackedDM(SupercellFns,SupercellFns),stat=ierr)
    if(ierr.ne.0) call stop_all(t_r,'Allocation error')
    UnpackedDM(:,:)=0.0_dp

    !Obtain unpacked DM

    !Just to test, use diagonal UnpackedDM
    UnpackedDM(:,:) = 0.0_dp
    do i=1,SupercellFns
        UnpackedDM(i,i) = 1.0_dp
    enddo
    
    !Extract the basis set over the supercell
    call c_f_pointer(Supercell%OrbBasis%Groups%pData, Supercell_groups, [Supercell%OrbBasis%Groups%nSize])
    !Extract the basis set over the unit cell 
    call c_f_pointer(Unitcell%OrbBasis%Groups%pData, Unitcell_groups, [Unitcell%OrbBasis%Groups%nSize])

    ngroups_super = Supercell%OrbBasis%Groups%nSize
    !write(6,*) "Number of groups of basis functions in supercell: ",ngroups_super
    ngroups_unit = Unitcell%OrbBasis%Groups%nSize
    !write(6,*) "Number of groups of basis functions in unitcell: ",ngroups_unit
    
    !write(6,*) "Unitcell lattice is: "
    !write(6,*) "Vector 1: ",Lattice%T(:,1)
    !write(6,*) "Vector 2: ",Lattice%T(:,2)
    !write(6,*) "Vector 3: ",Lattice%T(:,3)
!
!    write(6,*) "Supercell lattice is: "
!    write(6,*) "Vector 1: ",Supercell%T(:,1)
!    write(6,*) "Vector 2: ",Supercell%T(:,2)
!    write(6,*) "Vector 3: ",Supercell%T(:,3)
    
    !TODO: This should be called and set up once at the beginning
    call ConstructLattSumOrbPairs(Ex_Codens,Lattice,UnitCell,Supercell,CS_CodensThresh)

!In future, there should be a routine here for constructing the optimal Rc 
    !Want Rc .le. L/2
    !Calculate average L for supercell
    L_x = sqrt(Supercell%T(1,1)**2 + Supercell%T(2,1)**2 + Supercell%T(3,1)**2)
    L_y = sqrt(Supercell%T(1,2)**2 + Supercell%T(2,2)**2 + Supercell%T(3,2)**2)
    L_z = sqrt(Supercell%T(1,3)**2 + Supercell%T(2,3)**2 + Supercell%T(3,3)**2)
    L = (L_x + L_y + L_z) / 3.0_dp
    write(6,*) "Average length of supercell dimensions is: ",L
    Rc = L/2.0_dp 
    write(6,*) "Truncated coulomb potential kernel cutoff: ",Rc
    zTrans(:) = 0.0_dp  !The null translation vector

    ic = create_integral_context(0,0,1.0e-10_dp)
    !Plain coulomb
    !call assign_integral_kernel(ic,3,0,0)  !3 is for coulomb kernal
    !Truncated coulomb
    call assign_integral_kernel(ic,6,0,Rc)  !6 is for truncated coulomb kernal
    
    !What is the largest number of integrals (ij|ij) that can come out of two shells?
    !This is so we can store them when we calculate the exchange integrals.
    nlargestshell = 0
    do i = 1,ngroups_super
        if(Supercell_groups(i)%nFn.gt.nlargestshell) nlargestshell = Supercell_groups(i)%nFn
    enddo
    write(6,*) "Largest number of functions in a group: ",nlargestshell
    !How many integrals does this mean that we can get out of an (ij|ij) set?
    allocate(int_temp(nlargestshell,nlargestshell,nlargestshell,nlargestshell)) !This is to temporarily hold the integrals
    strides = (/ 1 , nlargestshell , nlargestshell**2 , nlargestshell**3 /)

    !Array Point to starting function of each shell - 1
    allocate(FunctionInd_Shell(ngroups_super))
    FunctionInd_Shell(:) = 0
    do i=2,ngroups_super
        FunctionInd_Shell(i) = FunctionInd_Shell(i-1) + Supercell_groups(i-1)%nFn
    enddo
!    write(6,*) "FunctionInd_Shell = ",FunctionInd_Shell(:)

!TODO: Sort out a proper bound for this - perhaps by looking at the largest estimated range between two functions & Rc
    iMaxCodensTrans = 6
    write(6,*) "Maximum number of translations in each direction: ",iMaxCodensTrans

    !Now, construct exchange matrix
    do mu_nu = 1, Ex_Codens%ngroup_pairs_UC
        !Run over codensities, with the first index in the unit cell, and second in supercell

        !Index of the shells mu and nu
        mu_ind = Ex_Codens%grouppairlist%group_labels(1,mu_nu)
        nu_ind = Ex_Codens%grouppairlist%group_labels(2,mu_nu)

        nbf_mu = Supercell_groups(mu_ind)%nFn
        nbf_nu = Supercell_groups(nu_ind)%nFn

        mu_nu_trans = Ex_Codens%grouppairlist%transVec(:,mu_nu)
        mu_nu_CS = Ex_Codens%grouppairlist%CS_Fac(mu_nu)

        do lam_sig = 1, Ex_Codens%ngroup_pairs
            
            !Now, *near-field* screening based on CS inequality.
            !No dependence on b translation
            !write(6,*) mu_nu,Ex_Codens%ngroup_pairs_UC,lam_sig,Ex_Codens%ngroup_pairs, &
            !    mu_nu_CS*Ex_Codens%grouppairlist%CS_Fac(lam_sig)
            if((mu_nu_CS*Ex_Codens%grouppairlist%CS_Fac(lam_sig)).lt.NF_ScreenTol) cycle

            lam_ind = Ex_Codens%grouppairlist%group_labels(1,lam_sig)
            sig_ind = Ex_Codens%grouppairlist%group_labels(2,lam_sig)

            !Ignore contribution if P^(nu,sig) is small
            tLattSum = .false.
            mu = FunctionInd_Shell(mu_ind)
            nu = FunctionInd_Shell(nu_ind)
            lam = FunctionInd_Shell(lam_ind)
            sig = FunctionInd_Shell(sig_ind)
            do nubf = nu+1,nu+nbf_nu
                do sigbf = sig+1,sig+nbf_sig
                    if(abs(UnpackedDM(nubf,sigbf)).gt.DM_ScreenTol) then
                        tLattSum = .true.
                        exit
                    endif
                enddo
            enddo
            if(.not.tLattSum) then
                !write(6,*) "Not performing sum for : ",nu+1,sig+1
                cycle
            endif

            nbf_lam = Supercell_groups(lam_ind)%nFn
            nbf_sig = Supercell_groups(sig_ind)%nFn

            lam_sig_trans = Ex_Codens%grouppairlist%transVec(:,lam_sig)

            !Run over all supercell translations 'b'
            do lxd = 1,2    !Which direction to go in
                do lx = 0,iMaxCodensTrans
                    if((lx.eq.0).and.(lxd.eq.2)) cycle
                    if(lxd.eq.1) then
                        x = lx
                    else
                        x = -lx
                    endif
                    do lyd = 1,2
                        do ly = 0,iMaxCodensTrans
                            if(lyd.eq.1) then
                                y = ly
                            else
                                y = -ly
                            endif
                            if((ly.eq.0).and.(lyd.eq.2)) cycle
                            do lzd = 1,2
                                do lz = 0,iMaxCodensTrans
                                    if((lz.eq.0).and.(lzd.eq.2)) cycle
                                    if(lzd.eq.1) then
                                        z = lz
                                    else
                                        z = -lz
                                    endif

                                    if((.not.tIncTrans).and.(z.ne.0).and.(z.ne.-1)) then
                                        !We didn't want any of the group pairs in the last translation,
                                        !And we have only gone out further now. Cycle through the rest of this lz set
                                        !This shouldn't result in any less pairs added to the list, and just be
                                        !and efficiency improvement
                                        cycle
                                        !TODO: We could also save how far we went out here, so that we didn't need to go out
                                        !further for the -lz direction? This should also be done for lx and ly.
                                        !Test this
                                        !Change to a do while loop moving from the outside
                                        !do while(abs(z).le.iMaxCodensTrans)
                                    endif

                                    !Calculate translation vector for groups b
                                    TransVec(:) = 0.0_dp    !3-vector of displacements for group b
                                    TransVec(:) = x*Supercell%T(:,1) + y*Supercell%T(:,2) + z*Supercell%T(:,3)
                                    BplusCTrans(:) = TransVec(:) + lam_sig_trans(:)
                                    tIncTrans = .false. !flag to indicate whether any contributions to this translation entered the exch
                                    !write(6,"(A,3I5)") "Considering SC translation: ",x,y,z

                                    !TODO: Screen here based on codensity distance
                                    !Can we do it earlier, before we even choose to include a particular lattice translation?

                                    int_temp(:,:,:,:) = 0.0_dp
                                    !Calculate set of integrals between all functions between these groups
                                    call eval_group_int2e_tra_incr(int_temp(1,1,1,1),strides(1),    &
                                        1.0_dp,mu_ind-1,zTrans(1),nu_ind-1,mu_nu_trans(1),lam_ind-1,TransVec(1), &
                                        sig_ind-1,BplusCTrans(1),Supercell%OrbBasis,ic)

                                    !Index of the specific functions given by the shell
                                    mu = FunctionInd_Shell(mu_ind)
                                    nu = FunctionInd_Shell(nu_ind)
                                    lam = FunctionInd_Shell(lam_ind)
                                    sig = FunctionInd_Shell(sig_ind)
            
                                    do mubf = 1,nbf_mu
                                        mu = mu + 1
                                        do nubf = 1,nbf_nu
                                            nu = nu + 1
                                            do lambf = 1,nbf_lam
                                                lam = lam + 1
                                                do sigbf = 1,nbf_sig
                                                    sig = sig + 1

                                                    if(.not.tIncTrans) then
                                                        if(int_temp(mubf,nubf,lambf,sigbf).gt.NF_ScreenTol) then
                                                            tIncTrans=.true.
                                                            !TODO: We need to include things symmetrically - we need to truncate BOTH
                                                            !sides equally....
                                                        endif
                                                    endif

                                                    write(66,*) x**2 + y**2 + z**2, int_temp(mubf,nubf,lambf,sigbf)

                                                    ExMat(mu,lam) = ExMat(mu,lam) + UnpackedDM(nu,sig)*   &
                                                        int_temp(mubf,nubf,lambf,sigbf)

                                                enddo
                                            enddo
                                        enddo
                                    enddo 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


!    nSS_Sq = SupercellFns*SupercellFns
!
!    !The number of functions in each shell is given by Supercell_groups(igroup)%nFn
!    mu=1
!    lam=1
!
!    !Loop over mu, lambda   SHELLS!
!    do mu_shell=1,ngroups_unit 
!        !Only want shells on symmetry unique atoms here for point group symmetry
!        nmu = Unitcell_groups(mu_shell)%nFn
!        write(6,*) 'unit cell fock index orbital shell: ',mu_shell,' / ',ngroups_unit
!
!        do lam_shell=1,ngroups_super
!            nlam = Supercell_groups(lam_shell)%nFn
!            write(6,*) 'supercell fock index orbital shell: ',lam_shell,' / ',ngroups_super
!
!            !Allocate memory for integrals converged over lattice vectors
!            allocate(ConvergedInts(SupercellFns,SupercellFns,nmu,nlam),stat=ierr)
!            if(ierr.ne.0) call stop_all(t_r,'Allocation error')
!            ConvergedInts(:,:,:,:) = 0.0_dp !Index 1: nu, index 2: sigma, index 3: mu (just over current group), 4: lambda (current group)
!            !Begin to accumulate (mu nu | lambda sigma) - nu fast
!
!            !up to 4 supercell translations in each direction
!            !This is specifically only for cubic symmetry - this needs to be generalized
!            !In non cubic symmetry, the shells are no longer 'degenerate', so perhaps we should just
!            !loop over in each direction. Probably best.
!            do shella=0,4 !16 
!                call findvectorsforshell(nperms_a,a_vecs,shella,48)
!                !write(6,*) "Looping over a shell: ",shella," Permutations: ",nperms_a
!                if(nperms_a.gt.48) call stop_all(t_r,'error in perms')
!                do ap=1,nperms_a
!
!                    do shellc=0,4 !16
!                        call findvectorsforshell(nperms_c,c_vecs,shellc,48)
!                        !write(6,*) "Looping over c shell: ",shellc," Permutations: ",nperms_c
!                        if(nperms_c.gt.48) call stop_all(t_r,'error in perms')
!                        do cp=1,nperms_c 
!
!                            do shellb=0,4 !16
!                                call findvectorsforshell(nperms_b,b_vecs,shellb,48)
!                                !write(6,*) "Looping over b shell: ",shellb," Permutations: ",nperms_b
!                                if(nperms_b.gt.48) call stop_all(t_r,'error in perms')
!                                do bp=1,nperms_b
!                                    !write(6,*) "b translation: ",bp,b_vecs(:,bp)
!
!                                    !The real-space translation vectors for the centers of the 
!                                    TransVec(:,1) = 0.0_dp  !Translation of mu (always in unit cell, therefore 0)
!                                    call GetSSTranslation(TransVec(:,2),a_vecs(:,ap),Supercell%T)      !Translation of nu (given by a_vecs(:,ap)
!                                    call GetSSTranslation(TransVec(:,3),b_vecs(:,bp),Supercell%T)      !Translation of lambda (given by b_vecs(:,bp)
!                                    call GetSSTranslation(TransVec(:,4),b_vecs(:,bp)+c_vecs(:,cp),Supercell%T) !Translation of sigma (b_vecs(:,bp)+c_vecs(:,cp))
!
!                                    !Now construct the integrals for all nu and sigma functions in the supercell
!                                    !These are not periodic integrals and the last three indices have their centers translated by various supercell vectors
!                                    !We should also be obtaining the integrals for the truncated coulomb potential
!                                    !(mu, nu [T a_vecs(:,ap)] | lambda [T b_vecs(:,bp)], sigma [(T b_vecs(:,bp)+c_vecs(:,cp)] ) 
!
!                                    ! Need to loop over all shells(/groups) of basis functions in supercell to construct matrix
!                                    !The basis functions are in supercell%OrbBasis
!
!                                    !Loop over nu, sigma   SHELLS!
!                                    sig=1
!                                    do sig_shell=1,ngroups_super
!                                        nsig = Supercell_groups(sig_shell)%nFn
!
!                                        nu=1
!                                        do nu_shell=1,ngroups_super
!                                            nnu = Supercell_groups(nu_shell)%nFn
!
!                                            !Want to fill up ConvergedInts(nu:nu+nnu,sig:sig+nsig,1:nmu,1:nlam) for (mu nu | lam sig)
!                                            strides = (/ nSS_sq, 1, nSS_sq*nMu, SupercellFns/) ! nSS x nSS x nMu x nLam
!                                            !write(6,*) "Loop over nu, sigma: ",nu,sig,ConvergedInts(nu,sig,1,1),sig_shell,lam_shell
!                                            !write(6,*) "TransVec: ",TransVec(:,:)
!                                            call eval_group_int2e_tra_incr(ConvergedInts(nu,sig,1,1),   &
!                                                        strides(1),-0.5_dp,mu_shell-1,TransVec(1,1),    &
!                                                        nu_shell-1,TransVec(1,2),lam_shell-1,TransVec(1,3), &
!                                                        sig_shell-1,TransVec(1,4),Supercell%OrbBasis, ic)
!
!                                            !write(6,*) "Loop over nu, sigma: ",nu,sig
!                                            !call flush(6)
!                                            !do i=1,SupercellFns
!                                            !    do j=1,SupercellFns
!                                            !        if(ConvergedInts(i,j,1,1).ne.0.0_dp) write(6,*) i,j,ConvergedInts(i,j,1,1)
!                                            !    enddo
!                                            !enddo
!                                            !write(6,*) ConvergedInts(:,:,:,:)
!                                            !call stop_all(t_r,"end")
!
!                                            nu=nu+nnu   !Increment the nu value we are up to
!                                        enddo
!                                        sig=sig+nsig    !Increment the sigma value we are up to
!                                    enddo
!
!
!                                enddo
!                            enddo
!                        enddo
!                    enddo
!                enddo
!            enddo
!
!
!            do xmu=1,nmu    !Loop over functions in mu shell
!                do xlam=1,nlam  !Loop over functions in lambda shell
!
!                    !contract here, to contract converged integrals with density matrix
!                    ExMat(mu+xmu-1,lam+xlam-1) = DDOT(nSS_sq,ConvergedInts(:,:,xmu,xlam),1,UnpackedDM,1)
!                    write(6,*) mu+xmu-1,lam+xlam-1,ExMat(mu+xmu-1,lam+xlam-1)
!
!                enddo
!            enddo
!
!            deallocate(ConvergedInts)
!
!            lam=lam+nlam
!        enddo
!        mu=mu+nmu
!    enddo
!
!    !Now calculate the exchange energy - contract again with the density matrix for the other electron
!    ExEnergy = 0.0_dp
!    ExEnergy = DDOT(nSS_Sq,ExMat,1,UnpackedDM,1)

    write(6,*) "Exchange matrix calculated."

end subroutine ExchangeSum

!!Get a set of translation vectors given by {dispvec(1)*SSvec(1), dispvec(2)*SSvec(2), dispvec(3)*SSvec(3)}
!pure subroutine GetSSTranslation(trans,dispvec,SSvec) 
!    use FortCons, only: dp
!    implicit none
!    integer, intent(in) :: dispvec(3)
!    real(dp), intent(in) :: SSvec(3,3)
!    real(dp), intent(out) :: trans(3)
!
!    trans(:) = dispvec(1)*SSvec(:,1) + dispvec(2)*SSvec(:,2) + dispvec(3)*SSvec(:,3)
!end subroutine GetSSTranslation
!
!!Return in vecs, the set of 'numvecs' set of {i,j,k} tuples, such that i^2+j^2+k^2=shell
!!This is written very specifically here - acually calculate it! Require recursion?
!subroutine findvectorsforshell(numvecs,vecs,shell,maxvecs)
!    use error_mod, only: stop_all
!    implicit none
!    integer, intent(in) :: shell,maxvecs
!    integer, intent(out) :: numvecs
!    integer, intent(out) :: vecs(3,maxvecs)
!    integer :: tuple(3),signedtuple(3)
!    integer :: sameint,sameint_signed,val
!    character(len=*), parameter :: t_r='findvectorsforshell'
!
!    if(shell.gt.16) call stop_all(t_r,'cannot deal with shells of distance > 4')
!
!    vecs(:,:)=0
!    tuple(:)=0
!    sameint=1
!    !Always order so that common pairs go last
!    if(shell.eq.0) then
!        sameint=3
!        continue
!    elseif(shell.eq.1) then
!        tuple=(/1,0,0/)
!        sameint=2
!    elseif(shell.eq.2) then
!        tuple=(/0,1,1/)
!        sameint=2
!    elseif(shell.eq.3) then
!        tuple=(/1,1,1/)
!        sameint=3
!    elseif(shell.eq.4) then
!        tuple=(/2,0,0/)
!        sameint=2
!    elseif(shell.eq.5) then
!        tuple=(/2,1,0/)
!    elseif(shell.eq.6) then
!        tuple=(/2,1,1/)
!        sameint=2
!    elseif(shell.eq.7) then
!        numvecs=0
!        return
!    elseif(shell.eq.8) then
!        tuple=(/0,2,2/)
!        sameint=2
!    elseif(shell.eq.9) then
!        tuple=(/3,0,0/)
!        sameint=2
!    elseif(shell.eq.10) then
!        tuple=(/3,1,0/)
!    elseif(shell.eq.11) then
!        tuple=(/3,1,1/)
!        sameint=2
!    elseif(shell.eq.12) then
!        tuple=(/2,2,2/) 
!        sameint=3
!    elseif(shell.eq.13) then
!        tuple=(/3,2,0/)
!    elseif(shell.eq.14) then
!        tuple=(/3,2,1/)
!    elseif(shell.eq.15) then
!        numvecs=0
!        return
!    elseif(shell.eq.16) then
!        tuple=(/4,0,0/)
!        sameint=2
!    endif
!
!    !Now, create all permuataions, including signed permutations
!    !Care with the fact that this may change sameint
!    val=1
!    ! + + +
!    call storeallperms(vecs,tuple,val,sameint,maxvecs)
!!    write(6,*) "val: ",val
!    if(tuple(1).ne.0) then
!        ! - + +
!        signedtuple(:)=tuple(:)
!        signedtuple(1)=-signedtuple(1)
!        if(sameint.ne.1) then
!            !i.e. some of the numbers are the same, therefore the number of permutations may be different
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 1')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "2 val: ",val
!    endif
!    if((tuple(2).ne.0).and.(tuple(1).ne.tuple(2))) then
!        ! + - +     !If the second value is the same as the first, we will have already counted these permutations
!        signedtuple(:)=tuple(:)
!        signedtuple(2)=-signedtuple(2)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) then
!                write(6,*) sameint_signed, sameint, signedtuple(:)
!                call stop_all(t_r,'error in calc perms 2')
!            endif
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "3 val: ",val
!    endif
!    if((tuple(3).ne.0).and.(tuple(3).ne.tuple(1)).and.(tuple(3).ne.tuple(2))) then
!        ! + + -
!        signedtuple(:)=tuple(:)
!        signedtuple(3)=-signedtuple(3)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 3')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "4 val: ",val
!    endif
!    if((tuple(1).ne.0).and.(tuple(2).ne.0)) then
!        ! - - +
!        signedtuple(:)=tuple(:)
!        signedtuple(1)=-signedtuple(1)
!        signedtuple(2)=-signedtuple(2)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 4')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "5 val: ",val
!    endif
!    if((tuple(1).ne.0).and.(tuple(3).ne.0).and.(sameint.eq.1)) then
!        !if all elements the same - ignore (sameint=3)
!        !if two elements the same - ignore, since always last two which are the same
!        ! - + -
!        signedtuple(:)=tuple(:)
!        signedtuple(1)=-signedtuple(1)
!        signedtuple(3)=-signedtuple(3)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 5')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "6 val: ",val
!    endif
!    if((tuple(2).ne.0).and.(tuple(3).ne.0).and.(sameint.ne.3)) then
!        ! + - -
!        signedtuple(:)=tuple(:)
!        signedtuple(2)=-signedtuple(2)
!        signedtuple(3)=-signedtuple(3)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 6')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "7 val: ",val
!    endif
!    if((tuple(1).ne.0).and.(tuple(2).ne.0).and.(tuple(3).ne.0)) then
!        ! - - -
!        signedtuple(:)=tuple(:)
!        signedtuple(1)=-signedtuple(1)
!        signedtuple(2)=-signedtuple(2)
!        signedtuple(3)=-signedtuple(3)
!        if(sameint.ne.1) then
!            call calcnumberpairs(sameint_signed,signedtuple)
!            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 7')
!        else
!            sameint_signed=sameint
!        endif
!        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!!        write(6,*) "8 val: ",val
!    endif
!
!    numvecs=val-1
!end subroutine findvectorsforshell
!
!pure subroutine calcnumberpairs(sameint,tuple)
!    implicit none
!    integer, intent(out) :: sameint
!    integer, intent(in) :: tuple(3)
!
!    if((tuple(1).eq.tuple(2)).and.(tuple(1).eq.tuple(3))) then
!        sameint=3   !all same
!    elseif((tuple(1).ne.tuple(2)).and.(tuple(1).ne.tuple(3)).and.(tuple(2).ne.tuple(3))) then
!        sameint=1   !all different
!    else
!        sameint=2   !one pair
!    endif
!end subroutine calcnumberpairs
!
!subroutine storeallperms(vecs,tuple,val,sameint,maxvecs)
!    implicit none
!    integer, intent(in) :: maxvecs
!    integer, intent(in) :: tuple(3) !the primitive permutation
!    integer, intent(in) :: sameint  !the number of pairs of identical integers
!    integer, intent(inout) :: val   !next free slot in vecs
!    integer, intent(inout) :: vecs(3,maxvecs)
!    
!    if(sameint.eq.1) then
!        !all integers unique
!        vecs(:,val)=tuple(:)                    !abc
!        call swap(vecs(:,val+1),vecs(:,val),2,3)     !acb
!        call swap(vecs(:,val+2),vecs(:,val+1),3,1)   !bca
!        call swap(vecs(:,val+3),vecs(:,val+2),2,3)   !bac
!        call swap(vecs(:,val+4),vecs(:,val+3),1,3)   !cab
!        call swap(vecs(:,val+5),vecs(:,val+4),2,3)   !cba
!        val=val+6
!    elseif(sameint.eq.2) then
!        !last two integers the same
!        vecs(:,val)=tuple(:)                    !abc
!        call swap(vecs(:,val+1),tuple(:),1,2)        !bac 
!        call swap(vecs(:,val+2),tuple(:),1,3)        !cba
!        val=val+3
!    elseif(sameint.eq.3) then
!        vecs(:,val)=tuple(:)
!        val=val+1
!    endif
!
!end subroutine storeallperms
!
!!Return the tuple outvec which is invec with the indices i and j permuted
!subroutine swap(outvec,invec,i,j)
!    implicit none
!    integer, intent(in) :: invec(3)
!    integer, intent(in) :: i,j
!    integer, intent(out) :: outvec(3)
!    integer :: k
!
!    outvec(:) = invec(:)
!    k=outvec(j)
!    outvec(j)=outvec(i)
!    outvec(i)=k
!end subroutine swap

