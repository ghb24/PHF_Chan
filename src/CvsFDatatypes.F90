module CvsFDatatypes
    use iso_c_binding
    implicit none

    integer, parameter :: sp = selected_real_kind(6,37)
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: qp = selected_real_kind(33,4931)
    integer, parameter :: int32 = selected_int_kind(8)
    integer, parameter :: int64 = selected_int_kind(15)

    type :: dyn_array_t
        type(c_ptr) :: pData
        integer :: nSize
        integer :: nReserved ! <-- don't touch that.
    end type

    type :: GroupInfo
        integer :: nFn
        integer :: nExp
        integer :: iExp
        integer :: nCo
        integer :: iCo
        integer :: l
        integer :: Type
        integer :: iCen
        integer :: iRange
    end type

    type :: basis_set_t
        type(dyn_array_t) :: Groups ! This is of type GroupInfo when unpacked
        type(dyn_array_t) :: Data   ! Type real
        type(dyn_array_t) :: Centers   !Type real(3)
        integer :: nFn
    end type

    type :: unit_cell_t
        type(dyn_array_t) :: Coords     !Type real(3)
        type(dyn_array_t) :: Elements   !Type Integer
        type(dyn_array_t) :: EcpCharges !Type Integer
        type(basis_set_t) :: OrbBasis
        real(dp) :: Volume
    end type

    type :: lattice_t
        !Real-space lattice vectors (translations between two unit-cells)
        real(dp), dimension(3,3) :: T
        !k-space lattice vectors (2pi*bi-orth basis of UnitCell(T)
        real(dp) , dimension(3,3) :: K
        real(dp) :: UnitCellVolume
    end type

    type :: supercell_t
        ! real-space displacements between two super-cells in the
        ! directions Lattice.T[i].
        ! This is the same as the translation between a raw Gaussian
        ! basis function and its periodic images.
        ! (note: SuperCell.T[i] == Size[i] * Lattice.T[i])
        real(dp), dimension(3,3) :: T
        ! [i]: number of times the unit-cell is repeated in Lattice.T[i]
        ! direction to form a super-cell. In total we obtain a
        !  Size[0] x Size[1] x Size[2]
        ! grid of super-cells.
        ! Note: Total number of super-cells is SuperCell.Ts.size().
        integer , dimension(3) :: Size
        ! set of all translation vectors between the first unit-cell
        ! and the other unit-cells in the super-cell (i.e., T x Size unpacked)
        ! Note: The 0-vector (first to first) is included in the set.
        type(dyn_array_t) :: Ts !Type real(3)
        type(basis_set_t) :: OrbBasis   
    end type

    type :: SolidModel
        type(unit_cell_t) :: UnitCell
        type(lattice_t) :: Lattice
        type(supercell_t) :: Supercell
    end type

end module
