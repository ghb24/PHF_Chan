module dummy_module
   use iso_c_binding
   type :: ints_and_floats_t
      double precision :: a, b
      integer :: x
      double precision :: c
      integer :: y
   end type

   type :: dyn_array_t
      type(c_ptr) :: pData
      integer :: nSize
      integer :: nReserved  ! <- don't touch that.
   end type

   type :: 

   type :: unit_cell_t
      type(dyn_array_t) :: Coords
      type(dyn_array_t) :: Elements
      type(dyn_array_t) :: EcpCharges
   end type

   contains

   subroutine aa(datastruc)
      type(unit_cell_t) :: datastruc
      integer, pointer :: Elements(:)
      call c_f_pointer(datastruc%Elements%pData, Elements, [datastruc%Elements%nSize])

end module

subroutine dummy_func_1(obj1)
   use dummy_module
   type(ints_and_floats_t) :: obj1
   obj1%x = -obj1%x
   obj1%c = -obj1%c
end subroutine

subroutine dummy_func_2(obj2)
   use dummy_module
   type(dyn_array_t) :: obj2
   integer, pointer :: Elements(:)
   call c_f_pointer(obj2%pData, Elements, [obj2%nSize])
   ! ^- these two tell fortran what kind of pointer we are supplying.
   !   note that since we can only declare pData as 'type(c_ptr)', the type
   !   information still present in c++ on what the arrays actually are are lost
   !   here and need to be given explicitly.
   !   I'm afraid I don't see a way around that, short of making allocatable arrays
   !   and TArray<> on C++ side binarily compatible... (would be very fragile)

   print '(1x,A,I4,A,I4)',  "(F90 side) size of the array: ", obj2%nSize, " reserved: ", obj2%nReserved
   print '(1x,A,99(1x,I10))', "(F90 side) elements:", (Elements(i),i=1,obj2%nSize)
end subroutine

subroutine dummy_func_3(obj2)
   use dummy_module
   type(dyn_array_t) :: obj2
   double precision, pointer :: Elements(:)
   call c_f_pointer(obj2%pData, Elements, [obj2%nSize])
   print '(1x,A,I4,A,I4)',  "(F90 side) size of the array: ", obj2%nSize, " reserved: ", obj2%nReserved
   print '(1x,A,99(1x,F10.4))', "(F90 side) elements:", (Elements(i),i=1,obj2%nSize)
end subroutine

