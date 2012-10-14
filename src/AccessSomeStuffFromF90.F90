

subroutine print_unit_cell(unit_cell)
   use cvsfdatatypes
   implicit double precision (a-h,o-z)
   type(unit_cell_t) :: unit_cell
   ! those are Fortran pointer refering to the actual pointers which are supplied in unit_cell
   double precision, pointer :: coords(:,:)
   integer, pointer :: elements(:)
   integer, pointer :: ecp_charges(:)
   ! we need to explicitly convert the pointers to tell Fortran what they are
   call c_f_pointer(unit_cell%coords%pdata, coords, [3,unit_cell%coords%nsize])
   call c_f_pointer(unit_cell%elements%pdata, elements, [unit_cell%elements%nsize])
   call c_f_pointer(unit_cell%ecpcharges%pdata, ecp_charges, [unit_cell%ecpcharges%nsize])

      print '(1x,A3,3x,A3,3x,3(2x,A10),4x,A3)',&
         '#AT', 'CHG', ' X ', ' Y ', ' Z ', 'ECP'
   do iAt = 1, size(elements)
      print '(1x,I3,3x,I3,3x,3(2x,F10.4),4x,I3)',&
         iAt, Elements(iAt), Coords(1,iAt), Coords(2,iAt), Coords(3,iAt), ecp_charges(iAt)
   end do
end subroutine