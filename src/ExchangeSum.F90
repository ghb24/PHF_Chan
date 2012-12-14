!Routine for constructing the a/c lattice sum list, and the orbital pairs that need to go with the
!exchange interaction.
!IN: Cell info
!    iMaxCodensTrans - maximum translation in any one direction from gamma-point to consider)
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

!Locally, also output largest translation vector in terms of sum of squared distance
subroutine ConstructLattSumOrbPairs(Ex_Codens,Lattice,UnitCell,Supercell,iMaxCodensTrans,dCoDensTol)
    use CvsFDatatypes
    use error_mod, only: stop_all
    use FortCons, only: dp
    implicit none
    !Arguments
    type(exchange_densities_t), intent(out) :: Ex_Codens
    type(lattice_t) , intent(in) :: Lattice
    type(unit_cell_t), intent(in)  :: UnitCell
    type(supercell_t), intent(in) :: Supercell
    integer, intent(in) :: iMaxCodensTrans
    real(dp), intent(in) :: dCoDensTol
    !Unpacked arrays
    type(GroupInfo), pointer :: Supercell_groups(:)
    type(GroupInfo), pointer :: Unitcell_groups(:)
    !local arguments
    integer :: UnitCellFns,SupercellFns,abf,bbf
    integer :: ngroups_super,ngroups_unit,nTrans,iNumGroups
    integer :: ic,nlargestshell,strides(4),iMaxTrans
    integer :: i,a,b,nbf_a,nbf_b,lx,lxd,x,ly,lyd,y,z,lz,lzd
    real(dp) :: L,Rc,zTrans(3),TransVec(3),magint
    integer, allocatable :: MaxTransVecs(:,:)
    real(dp), allocatable :: int_temp(:,:,:,:)
    logical ::  tIncTrans
    !external functions
    integer :: create_integral_context
    
    UnitCellFns = UnitCell%OrbBasis%nFn
    SupercellFns = Supercell%OrbBasis%nFn
    
    !Extract the basis set over the supercell
    call c_f_pointer(Supercell%OrbBasis%Groups%pData, Supercell_groups, [Supercell%OrbBasis%Groups%nSize])
    !Extract the basis set over the unit cell 
    call c_f_pointer(Unitcell%OrbBasis%Groups%pData, Unitcell_groups, [Unitcell%OrbBasis%Groups%nSize])

    ngroups_super = Supercell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in supercell: ",ngroups_super
    write(6,*) "Number of functions in supercell: ",SupercellFns
    ngroups_unit = Unitcell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in unitcell: ",ngroups_unit
    
!In future, there should be a routine here for constructing the optimal Rc 
    !Want Rc .le. L/2
    !Calculate L for cubic cell - just take one direction. Won't be right in general
    L = sqrt(Supercell%T(1,1)**2 + Supercell%T(2,1)**2 + Supercell%T(3,1)**2)
    write(6,*) "Length of supercell is approximately: ",L
    Rc = L/2.0_dp 
    write(6,*) "Truncated coulomb potential kernel cutoff: ",Rc

    ic = create_integral_context(0,0,1.0e-10_dp)
    !Plain coulomb
    !call assign_integral_kernel(ic,3,0,0)  !3 is for coulomb kernal
    !Truncated coulomb
    call assign_integral_kernel(ic,6,0,Rc)  !6 is for truncated coulomb kernal

    write(6,"(A)") "Constructing list of codensities for exchange summation"

    !What is the largest number of integrals (ij|ij) that can come out of two shells?
    nlargestshell = 0
    do i = 1,ngroups_super
        if(Supercell_groups(i)%nFn.gt.nlargestshell) nlargestshell = Supercell_groups(i)%nFn
    enddo
    write(6,*) "Largest number of functions in a group: ",nlargestshell
    !How many integrals does this mean that we can get out of an (ij|ij) set?
    allocate(int_temp(nlargestshell,nlargestshell,nlargestshell,nlargestshell)) !This is to temporarily hold the integrals
    strides = (/ 1 , nlargestshell , nlargestshell**2 , nlargestshell**3 /)

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
    do lx = 0,iMaxCodensTrans
        do lxd = 1,2    !Which direction to go in
            if((lx.eq.0).and.(lxd.eq.2)) cycle
            if(lxd.eq.1) then
                x = lx
            else
                x = -lx
            endif
            do ly = 0,iMaxCodensTrans
                do lyd = 1,2
                    if(lyd.eq.1) then
                        y = ly
                    else
                        y = -ly
                    endif
                    if((ly.eq.0).and.(lyd.eq.2)) cycle
                    do lz = 0,iMaxCodensTrans
                        do lzd = 1,2
                            if((lz.eq.0).and.(lzd.eq.2)) cycle
                            if(lzd.eq.1) then
                                z = lz
                            else
                                z = -lz
                            endif

                            !Calculate translation vector for groups b
                            TransVec(:) = 0.0_dp    !3-vector of displacements for group b
                            TransVec(:) = x*Supercell%T(:,1) + y*Supercell%T(:,2) + z*Supercell%T(:,3)
                            tIncTrans = .false. !flag to indicate whether to include this translation in the list
                            write(6,"(A,3I5)") "Considering SC translation: ",x,y,z

                            do a = 1,ngroups_super

                                write(6,*) "Considering gamma-point group: ",a
                                nbf_a = Supercell_groups(a)%nFn

                                !Now loop over functions in this translated supercell
                                do b = 1,ngroups_super

                                    nbf_b = Supercell_groups(b)%nFn
                                    int_temp(:,:,:,:) = 0.0_dp
                                    !Calculate set of integrals between all functions between these groups
                                    call eval_group_int2e_tra_incr(int_temp(1,1,1,1),strides(1),    &
                                        1.0_dp,a-1,zTrans(1),b-1,TransVec(1),a-1,zTrans(1),b-1,     &
                                        TransVec(1),Supercell%OrbBasis,ic)

                                    !TODO: Check these orbitals are now packed as expected

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
                                                write(16,*) x**2+y**2+z**2,magint
                                                !exit loop
                                            endif
                                        enddo
                                    enddo loop

                                enddo   !end b

                                !TODO: We need an exit criterion to get out of these loops if we have definitely gone far enough

                            enddo   !lzd
                        enddo   !lz
                    enddo   !lyd
                enddo   !ly
            enddo   !lxd
        enddo   !lx
    enddo   ! end a

    write(6,*) "SC Translation vectors that we need to consider for codensities: "
    do i = 1,nTrans
        write(6,"(3I6)") MaxTransVecs(:,i)
    enddo
    write(6,"(A,I9,A,I9,A)") "This corresponds to ",nTrans," SC translations out of a possible searched ",iMaxTrans," translations."
    write(6,*) ""
    write(6,"(A,I15)") "Total number of group pairs to consider for each electron in exchange sum: ",iNumGroups

    !Now, move this list into the appropriate derived type
    !allocate(Ex_Codens%TransVecs%pData(3,nTrans))
    !Ex_Codens%TransVecs%nSize = nTrans*3

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

    !test
    type(Exchange_densities_t) :: Ex_Codens

    !Unpacked arguments
    real(dp), pointer :: DenMat(:,:)
    real(dp), pointer :: ExMat(:,:)
    type(GroupInfo), pointer :: Supercell_groups(:)
    type(GroupInfo), pointer :: Unitcell_groups(:)

    !Local variables
    real(dp) :: dDenDim   !Dimension of density matrix
    real(dp) :: Rc,L
    integer :: UnitCellFns
    integer :: SupercellFns
    integer :: ap,bp,cp,i,j,ic,ierr
    integer :: nnu,nsig,nlam,nmu,ngroups_unit,ngroups_super,nperms_a,nperms_b,nperms_c
    integer :: lam_shell,mu_shell,sig_shell,nu_shell
    integer :: lam,mu,nu,sig,xlam,xmu,nSS_Sq
    integer :: shella,shellb,shellc,SymUniqUnitcellFns
    real(dp), allocatable :: UnpackedDM(:,:)    !SS x SS unpacked DM
    real(dp), allocatable :: ConvergedInts(:,:,:,:) !The converged integrals
    !The lattice translations for the sums
    integer, allocatable :: a_vecs(:,:),b_vecs(:,:),c_vecs(:,:)   
    !Real-space translations of each orbital in the integral sum
    real(dp) :: TransVec(3,4)
    integer :: strides(4)
    character(len=*), parameter :: t_r='ExchangeSum'

    !external functions
    integer :: create_integral_context
    real(dp) :: DDOT

    call ConstructLattSumOrbPairs(Ex_Codens,Lattice,UnitCell,Supercell,6,1.0e-8_dp)

    call stop_all(t_r,'end')

    UnitCellFns = UnitCell%OrbBasis%nFn
    write(6,*) "Number of basis functions in unit cell: ",UnitCellFns
    SupercellFns = Supercell%OrbBasis%nFn
    write(6,*) "Number of basis functions in supercell: ",SupercellFns

    write(6,*) "Number of unit cells in supercell: ",Supercell%Size(:)
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

    write(6,*) "Size of matrices as expected"

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

    !There should be a maximum of 48 permutations for any one shell, assuming no degeneracy
    allocate(a_vecs(3,48))
    allocate(b_vecs(3,48))
    allocate(c_vecs(3,48))

    !Extract the basis set over the supercell
    call c_f_pointer(Supercell%OrbBasis%Groups%pData, Supercell_groups, [Supercell%OrbBasis%Groups%nSize])
    !Extract the basis set over the unit cell 
    call c_f_pointer(Unitcell%OrbBasis%Groups%pData, Unitcell_groups, [Unitcell%OrbBasis%Groups%nSize])

    ngroups_super = Supercell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in supercell: ",ngroups_super
    ngroups_unit = Unitcell%OrbBasis%Groups%nSize
    write(6,*) "Number of groups of basis functions in unitcell: ",ngroups_unit
    
    !write(6,*) "Unitcell lattice is: "
    !write(6,*) "Vector 1: ",Lattice%T(:,1)
    !write(6,*) "Vector 2: ",Lattice%T(:,2)
    !write(6,*) "Vector 3: ",Lattice%T(:,3)
!
!    write(6,*) "Supercell lattice is: "
!    write(6,*) "Vector 1: ",Supercell%T(:,1)
!    write(6,*) "Vector 2: ",Supercell%T(:,2)
!    write(6,*) "Vector 3: ",Supercell%T(:,3)

    !Want Rc .le. L/2
    !Calculate L for cubic cell - just take one direction. Won't be right in general
    L = sqrt(Supercell%T(1,1)**2 + Supercell%T(2,1)**2 + Supercell%T(3,1)**2)
    write(6,*) "Length of supercell is approximately: ",L
    Rc = L/2.0_dp 
    write(6,*) "Truncated coulomb potential kernel cutoff: ",Rc

    ic = create_integral_context(0,0,1.0e-10_dp)
    !Plain coulomb
    !call assign_integral_kernel(ic,3,0,0)  !3 is for coulomb kernal
    !Truncated coulomb
    call assign_integral_kernel(ic,6,0,Rc)  !6 is for truncated coulomb kernal

    nSS_Sq = SupercellFns*SupercellFns

    !The number of functions in each shell is given by Supercell_groups(igroup)%nFn
    mu=1
    lam=1

    !Loop over mu, lambda   SHELLS!
    do mu_shell=1,ngroups_unit 
        !Only want shells on symmetry unique atoms here for point group symmetry
        nmu = Unitcell_groups(mu_shell)%nFn
        write(6,*) 'unit cell fock index orbital shell: ',mu_shell,' / ',ngroups_unit

        do lam_shell=1,ngroups_super
            nlam = Supercell_groups(lam_shell)%nFn
            write(6,*) 'supercell fock index orbital shell: ',lam_shell,' / ',ngroups_super

            !Allocate memory for integrals converged over lattice vectors
            allocate(ConvergedInts(SupercellFns,SupercellFns,nmu,nlam),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Allocation error')
            ConvergedInts(:,:,:,:) = 0.0_dp !Index 1: nu, index 2: sigma, index 3: mu (just over current group), 4: lambda (current group)
            !Begin to accumulate (mu nu | lambda sigma) - nu fast

            !up to 4 supercell translations in each direction
            !This is specifically only for cubic symmetry - this needs to be generalized
            !In non cubic symmetry, the shells are no longer 'degenerate', so perhaps we should just
            !loop over in each direction. Probably best.
            do shella=0,4 !16 
                call findvectorsforshell(nperms_a,a_vecs,shella,48)
                !write(6,*) "Looping over a shell: ",shella," Permutations: ",nperms_a
                if(nperms_a.gt.48) call stop_all(t_r,'error in perms')
                do ap=1,nperms_a

                    do shellc=0,4 !16
                        call findvectorsforshell(nperms_c,c_vecs,shellc,48)
                        !write(6,*) "Looping over c shell: ",shellc," Permutations: ",nperms_c
                        if(nperms_c.gt.48) call stop_all(t_r,'error in perms')
                        do cp=1,nperms_c 

                            do shellb=0,4 !16
                                call findvectorsforshell(nperms_b,b_vecs,shellb,48)
                                !write(6,*) "Looping over b shell: ",shellb," Permutations: ",nperms_b
                                if(nperms_b.gt.48) call stop_all(t_r,'error in perms')
                                do bp=1,nperms_b
                                    !write(6,*) "b translation: ",bp,b_vecs(:,bp)

                                    !The real-space translation vectors for the centers of the 
                                    TransVec(:,1) = 0.0_dp  !Translation of mu (always in unit cell, therefore 0)
                                    call GetSSTranslation(TransVec(:,2),a_vecs(:,ap),Supercell%T)      !Translation of nu (given by a_vecs(:,ap)
                                    call GetSSTranslation(TransVec(:,3),b_vecs(:,bp),Supercell%T)      !Translation of lambda (given by b_vecs(:,bp)
                                    call GetSSTranslation(TransVec(:,4),b_vecs(:,bp)+c_vecs(:,cp),Supercell%T) !Translation of sigma (b_vecs(:,bp)+c_vecs(:,cp))

                                    !Now construct the integrals for all nu and sigma functions in the supercell
                                    !These are not periodic integrals and the last three indices have their centers translated by various supercell vectors
                                    !We should also be obtaining the integrals for the truncated coulomb potential
                                    !(mu, nu [T a_vecs(:,ap)] | lambda [T b_vecs(:,bp)], sigma [(T b_vecs(:,bp)+c_vecs(:,cp)] ) 

                                    ! Need to loop over all shells(/groups) of basis functions in supercell to construct matrix
                                    !The basis functions are in supercell%OrbBasis

                                    !Loop over nu, sigma   SHELLS!
                                    sig=1
                                    do sig_shell=1,ngroups_super
                                        nsig = Supercell_groups(sig_shell)%nFn

                                        nu=1
                                        do nu_shell=1,ngroups_super
                                            nnu = Supercell_groups(nu_shell)%nFn

                                            !Want to fill up ConvergedInts(nu:nu+nnu,sig:sig+nsig,1:nmu,1:nlam) for (mu nu | lam sig)
                                            strides = (/ nSS_sq, 1, nSS_sq*nMu, SupercellFns/) ! nSS x nSS x nMu x nLam
                                            !write(6,*) "Loop over nu, sigma: ",nu,sig,ConvergedInts(nu,sig,1,1),sig_shell,lam_shell
                                            !write(6,*) "TransVec: ",TransVec(:,:)
                                            call eval_group_int2e_tra_incr(ConvergedInts(nu,sig,1,1),   &
                                                        strides(1),-0.5_dp,mu_shell-1,TransVec(1,1),    &
                                                        nu_shell-1,TransVec(1,2),lam_shell-1,TransVec(1,3), &
                                                        sig_shell-1,TransVec(1,4),Supercell%OrbBasis, ic)

                                            !write(6,*) "Loop over nu, sigma: ",nu,sig
                                            !call flush(6)
                                            !do i=1,SupercellFns
                                            !    do j=1,SupercellFns
                                            !        if(ConvergedInts(i,j,1,1).ne.0.0_dp) write(6,*) i,j,ConvergedInts(i,j,1,1)
                                            !    enddo
                                            !enddo
                                            !write(6,*) ConvergedInts(:,:,:,:)
                                            !call stop_all(t_r,"end")

                                            nu=nu+nnu   !Increment the nu value we are up to
                                        enddo
                                        sig=sig+nsig    !Increment the sigma value we are up to
                                    enddo


                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo


            do xmu=1,nmu    !Loop over functions in mu shell
                do xlam=1,nlam  !Loop over functions in lambda shell

                    !contract here, to contract converged integrals with density matrix
                    ExMat(mu+xmu-1,lam+xlam-1) = DDOT(nSS_sq,ConvergedInts(:,:,xmu,xlam),1,UnpackedDM,1)
                    write(6,*) mu+xmu-1,lam+xlam-1,ExMat(mu+xmu-1,lam+xlam-1)

                enddo
            enddo

            deallocate(ConvergedInts)

            lam=lam+nlam
        enddo
        mu=mu+nmu
    enddo

    !Now calculate the exchange energy - contract again with the density matrix for the other electron
    ExEnergy = 0.0_dp
    ExEnergy = DDOT(nSS_Sq,ExMat,1,UnpackedDM,1)

    write(6,*) "Exchange energy calculated as: ",ExEnergy

end subroutine ExchangeSum

!Get a set of translation vectors given by {dispvec(1)*SSvec(1), dispvec(2)*SSvec(2), dispvec(3)*SSvec(3)}
pure subroutine GetSSTranslation(trans,dispvec,SSvec) 
    use FortCons, only: dp
    implicit none
    integer, intent(in) :: dispvec(3)
    real(dp), intent(in) :: SSvec(3,3)
    real(dp), intent(out) :: trans(3)

    trans(:) = dispvec(1)*SSvec(:,1) + dispvec(2)*SSvec(:,2) + dispvec(3)*SSvec(:,3)
end subroutine GetSSTranslation

!Return in vecs, the set of 'numvecs' set of {i,j,k} tuples, such that i^2+j^2+k^2=shell
!This is written very specifically here - acually calculate it! Require recursion?
subroutine findvectorsforshell(numvecs,vecs,shell,maxvecs)
    use error_mod, only: stop_all
    implicit none
    integer, intent(in) :: shell,maxvecs
    integer, intent(out) :: numvecs
    integer, intent(out) :: vecs(3,maxvecs)
    integer :: tuple(3),signedtuple(3)
    integer :: sameint,sameint_signed,val
    character(len=*), parameter :: t_r='findvectorsforshell'

    if(shell.gt.16) call stop_all(t_r,'cannot deal with shells of distance > 4')

    vecs(:,:)=0
    tuple(:)=0
    sameint=1
    !Always order so that common pairs go last
    if(shell.eq.0) then
        sameint=3
        continue
    elseif(shell.eq.1) then
        tuple=(/1,0,0/)
        sameint=2
    elseif(shell.eq.2) then
        tuple=(/0,1,1/)
        sameint=2
    elseif(shell.eq.3) then
        tuple=(/1,1,1/)
        sameint=3
    elseif(shell.eq.4) then
        tuple=(/2,0,0/)
        sameint=2
    elseif(shell.eq.5) then
        tuple=(/2,1,0/)
    elseif(shell.eq.6) then
        tuple=(/2,1,1/)
        sameint=2
    elseif(shell.eq.7) then
        numvecs=0
        return
    elseif(shell.eq.8) then
        tuple=(/0,2,2/)
        sameint=2
    elseif(shell.eq.9) then
        tuple=(/3,0,0/)
        sameint=2
    elseif(shell.eq.10) then
        tuple=(/3,1,0/)
    elseif(shell.eq.11) then
        tuple=(/3,1,1/)
        sameint=2
    elseif(shell.eq.12) then
        tuple=(/2,2,2/) 
        sameint=3
    elseif(shell.eq.13) then
        tuple=(/3,2,0/)
    elseif(shell.eq.14) then
        tuple=(/3,2,1/)
    elseif(shell.eq.15) then
        numvecs=0
        return
    elseif(shell.eq.16) then
        tuple=(/4,0,0/)
        sameint=2
    endif

    !Now, create all permuataions, including signed permutations
    !Care with the fact that this may change sameint
    val=1
    ! + + +
    call storeallperms(vecs,tuple,val,sameint,maxvecs)
!    write(6,*) "val: ",val
    if(tuple(1).ne.0) then
        ! - + +
        signedtuple(:)=tuple(:)
        signedtuple(1)=-signedtuple(1)
        if(sameint.ne.1) then
            !i.e. some of the numbers are the same, therefore the number of permutations may be different
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 1')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "2 val: ",val
    endif
    if((tuple(2).ne.0).and.(tuple(1).ne.tuple(2))) then
        ! + - +     !If the second value is the same as the first, we will have already counted these permutations
        signedtuple(:)=tuple(:)
        signedtuple(2)=-signedtuple(2)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) then
                write(6,*) sameint_signed, sameint, signedtuple(:)
                call stop_all(t_r,'error in calc perms 2')
            endif
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "3 val: ",val
    endif
    if((tuple(3).ne.0).and.(tuple(3).ne.tuple(1)).and.(tuple(3).ne.tuple(2))) then
        ! + + -
        signedtuple(:)=tuple(:)
        signedtuple(3)=-signedtuple(3)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 3')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "4 val: ",val
    endif
    if((tuple(1).ne.0).and.(tuple(2).ne.0)) then
        ! - - +
        signedtuple(:)=tuple(:)
        signedtuple(1)=-signedtuple(1)
        signedtuple(2)=-signedtuple(2)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 4')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "5 val: ",val
    endif
    if((tuple(1).ne.0).and.(tuple(3).ne.0).and.(sameint.eq.1)) then
        !if all elements the same - ignore (sameint=3)
        !if two elements the same - ignore, since always last two which are the same
        ! - + -
        signedtuple(:)=tuple(:)
        signedtuple(1)=-signedtuple(1)
        signedtuple(3)=-signedtuple(3)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 5')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "6 val: ",val
    endif
    if((tuple(2).ne.0).and.(tuple(3).ne.0).and.(sameint.ne.3)) then
        ! + - -
        signedtuple(:)=tuple(:)
        signedtuple(2)=-signedtuple(2)
        signedtuple(3)=-signedtuple(3)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 6')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "7 val: ",val
    endif
    if((tuple(1).ne.0).and.(tuple(2).ne.0).and.(tuple(3).ne.0)) then
        ! - - -
        signedtuple(:)=tuple(:)
        signedtuple(1)=-signedtuple(1)
        signedtuple(2)=-signedtuple(2)
        signedtuple(3)=-signedtuple(3)
        if(sameint.ne.1) then
            call calcnumberpairs(sameint_signed,signedtuple)
            if(sameint_signed.gt.sameint) call stop_all(t_r,'error in calc perms 7')
        else
            sameint_signed=sameint
        endif
        call storeallperms(vecs,signedtuple,val,sameint_signed,maxvecs)
!        write(6,*) "8 val: ",val
    endif

    numvecs=val-1
end subroutine findvectorsforshell

pure subroutine calcnumberpairs(sameint,tuple)
    implicit none
    integer, intent(out) :: sameint
    integer, intent(in) :: tuple(3)

    if((tuple(1).eq.tuple(2)).and.(tuple(1).eq.tuple(3))) then
        sameint=3   !all same
    elseif((tuple(1).ne.tuple(2)).and.(tuple(1).ne.tuple(3)).and.(tuple(2).ne.tuple(3))) then
        sameint=1   !all different
    else
        sameint=2   !one pair
    endif
end subroutine calcnumberpairs

subroutine storeallperms(vecs,tuple,val,sameint,maxvecs)
    implicit none
    integer, intent(in) :: maxvecs
    integer, intent(in) :: tuple(3) !the primitive permutation
    integer, intent(in) :: sameint  !the number of pairs of identical integers
    integer, intent(inout) :: val   !next free slot in vecs
    integer, intent(inout) :: vecs(3,maxvecs)
    
    if(sameint.eq.1) then
        !all integers unique
        vecs(:,val)=tuple(:)                    !abc
        call swap(vecs(:,val+1),vecs(:,val),2,3)     !acb
        call swap(vecs(:,val+2),vecs(:,val+1),3,1)   !bca
        call swap(vecs(:,val+3),vecs(:,val+2),2,3)   !bac
        call swap(vecs(:,val+4),vecs(:,val+3),1,3)   !cab
        call swap(vecs(:,val+5),vecs(:,val+4),2,3)   !cba
        val=val+6
    elseif(sameint.eq.2) then
        !last two integers the same
        vecs(:,val)=tuple(:)                    !abc
        call swap(vecs(:,val+1),tuple(:),1,2)        !bac 
        call swap(vecs(:,val+2),tuple(:),1,3)        !cba
        val=val+3
    elseif(sameint.eq.3) then
        vecs(:,val)=tuple(:)
        val=val+1
    endif

end subroutine storeallperms

!Return the tuple outvec which is invec with the indices i and j permuted
subroutine swap(outvec,invec,i,j)
    implicit none
    integer, intent(in) :: invec(3)
    integer, intent(in) :: i,j
    integer, intent(out) :: outvec(3)
    integer :: k

    outvec(:) = invec(:)
    k=outvec(j)
    outvec(j)=outvec(i)
    outvec(i)=k
end subroutine swap

