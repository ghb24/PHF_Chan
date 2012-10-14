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

    !Unpacked arguments
    real(dp), pointer :: DenMat(:,:)
    real(dp), pointer :: ExMat(:,:)

    !Local variables
    real(dp) :: dDenDim   !Dimension of density matrix
    integer :: UnitCellFns
    integer :: SupercellFns
    character(len=*), parameter :: t_r='ExchangeSum'

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

end subroutine ExchangeSum
