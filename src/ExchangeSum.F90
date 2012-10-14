subroutine ExchangeSum(Lattice,UnitCell)
    use CvsFDatatypes
    type(lattice_t) , intent(in) :: Lattice
    type(unit_cell_t), intent(in)  :: UnitCell

    integer, pointer :: Elements(:) !Extract the fortran data from the c pointer
                                    !Is this of dimension 2?

    call c_f_pointer(UnitCell%Elements%pData, Elements, [UnitCell%Elements%nSize])

    write(6,*) "Elements in cell are: ",Elements(:)

end subroutine ExchangeSum
