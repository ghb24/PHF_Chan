module CvsFDatatypes
    use iso_c_binding
    use FortCons, only: dp
    implicit none

    !This is the fortran array for passing dynamic arrays from C objects
    !pData contains a pointer to the actual array, where nSize is the size of the array
    type, bind(c) :: dyn_array_t
        type(c_ptr) :: pData
        integer :: nSize
        integer :: nReserved ! <-- don't touch that.
    end type

! A note on Gaussian basis sets (general):
!  - A single Gaussian basis basis function (e.g., a 2px function) cannot be
!    evaluated efficiently. Normally, there are many functions which share the
!    same evaluation intermediates (e.g., p functions on a center); these have
!    to be evaluated together
!
!  - For this reason a basis set is subdivided into /shell groups/. A shell
!    group bundles all functions using shared intermediates. Normally a shell
!    group consists of:
!
!      o Functions on a single center, with the same angular
!        momentum
!      o One set of primitive exponents used by all of them
!      o A contraction matrix defining how the primitives are
!        linearly combined to form the actual basis functions.
!
!    In rare cases (the one we're about to implement being one) it is sensible
!    to share exponents between the s- and p functions (and thus create
!    dual-angular momentum shell groups). For simplicitly reasons this
!    is not done here.
!
!  - A shell group has (2*l + 1)*nCo functions, where 'l' is the angular
!    momentum and nCo is the number of contractions.
!
! In the solid case:
!  - For solids, technically the basis functions *themselves* are periodic,
!    with the period of the super-cell size (see supercells_and_kspace.pdf).
!
!  - I(cgk)'ll try to make the integral/basis-function-on-grid routines etc
!    take that into account automatically. For this reason the basis set
!    objects carry super-cell periodicity vectors.
!
!
! note:
!  - the objects here and on the Fortran side need to be consistent!
!
!  - Within this program, all indices are 0-based!
    type, bind(c) :: GroupInfo
        integer :: nFn  !Number of functions (contractions x components)
        integer :: nExp !Number of primitive exponents
        integer :: iExp !Start index in exponent array in externally 
                        !supplied *pExp / *pBasisData
        integer :: nCo  !Number of contractions
        integer :: iCo  !Start index of nExp x nCo contraction matrix, 
                        !in externally supplied *pCo / *pBasisData
        integer :: l    !Angular momentum
        integer :: Type !Function type (0=GTO, 1=Poisson)
        integer :: iCen !Center index
        integer :: iRange   !pBasisData[i]: effective range at which we 
                            !can consider the basis function to be zero
    end type

!Holds information on the gaussian basis set
    type, bind(c) :: basis_set_t
        type(dyn_array_t) :: Groups ! This is of type GroupInfo when unpacked
                                    ! Shell groups: defines the actual basis functions
        type(dyn_array_t) :: Data   ! Type real
                                    ! Data on sequences of (flattened) contraction matrices 
                                    ! and exponent arrays. Data in *pGroupInfo refers to these.
        type(dyn_array_t) :: Centers   !Type real(3)
                                       ! Centers of the raw gaussians
        integer :: nFn              !Total number of basis functions
    end type

    type, bind(c) :: unit_cell_t
        type(dyn_array_t) :: Coords     !Type real(3), xyz coordinates of atoms
        type(dyn_array_t) :: Elements   !Type Integer, atom types
        type(dyn_array_t) :: EcpCharges !Type Integer, electrons absorbed into PP
        type(basis_set_t) :: OrbBasis   !Gaussian basis specification with centers in the unitcell
                                        !Actual basis is obv periodically repeated with the super-cell repetition
        real(dp) :: Volume              !3d volume (=det(Lattice.T))
    end type
    
    !Type for codensity of two group functions
    type :: group_pairs_t
        integer, allocatable :: group_labels(:,:)   !(2,ngroup_pairs) - List of group pairs to consider in lattice sums, labelled by supercell indices
        real(dp), allocatable :: transVec(:,:)       !(3,ngroup_pairs) - real translation vector of second group index
        real(dp), allocatable :: CS_Fac(:)          !(ngroup_pairs) - max value of |(ij_c|ij_c)|^0.5 between the functions in the two groups
        real(dp), allocatable :: den_radius(:)      !(ngroup_pairs) - approximate max radius of the codensity between functions in the two groups
        real(dp), allocatable :: den_center(:,:)    !(3,ngroup_pairs) - approximate center of the codensity between functions in the two groups
    end type

    !Type which corresponds to the list of codensities to consider in the exchange sum
    type :: exchange_densities_t
        integer :: nTransVecs                   !Number of unique translation vectors between all pairs of groups to consider
        integer, allocatable :: iTransVecs(:,:)  !(3,nTransVecs) - supercell translation factors to consider
        integer :: ngroup_pairs                 !Number of group pairs to consider in codensity sum
        integer :: ngroup_pairs_UC              !Number of group pairs where first group index is in the unit cell (these are indexed first)
        type(group_pairs_t) :: grouppairlist    !List of group pairs
    end type

! defines the physical lattice: allowed periodic translations of the unit-cell.
!(note: Lattice symmetry stuff, classical lattice summations(*), and lattice
! reduction algorithms should probably go here.
! (*): That is a routine which takes a number of point charges in the unit
!      cell (neutral in total), and a (large) number of real-space grid
!      points, and evaluates the classical Madelung potential on the grid.)
    type, bind(c) :: lattice_t
        !Real-space lattice vectors (translations between two unit-cells)
        real(dp), dimension(3,3) :: T
        !k-space lattice vectors (2pi*bi-orth basis of UnitCell(T)
        real(dp) , dimension(3,3) :: K
        real(dp) :: UnitCellVolume
    end type

! defines our calculation model: number of unit-cells explicitly treated
! before being periodically repeated.
! Notes:
!   - Unit-cell placement: First unit-cell is at T = (0,0,0), and other
!     unit-cells are displaced in *positive* Lattice.T direction (i.e.,
!     the first unit-cell is at the super-cell boundary).
!     This might simplify some implementation aspects, and theoretically
!     it should not matter since the calculation is periodic)
!   - As the k-point grid is effectively a different representation of
!     the super-cell information, it should probably also go in here
!     once we get to that point.
    type, bind(c) :: supercell_t
        real(dp), dimension(3,3) :: T
        ! real-space displacements between two super-cells in the
        ! directions Lattice.T[i].
        ! This is the same as the translation between a raw Gaussian
        ! basis function and its periodic images.
        ! (note: SuperCell.T[i] == Size[i] * Lattice.T[i])
        integer , dimension(3) :: Size
        ! [i]: number of times the unit-cell is repeated in Lattice.T[i]
        ! direction to form a super-cell. In total we obtain a
        !  Size[0] x Size[1] x Size[2]
        ! grid of super-cells.
        ! Note: Total number of super-cells is SuperCell.Ts.size().
        type(dyn_array_t) :: Ts !Type real(3)
        ! set of all translation vectors between the first unit-cell
        ! and the other unit-cells in the super-cell (i.e., T x Size unpacked)
        ! Note: The 0-vector (first to first) is included in the set.
        type(basis_set_t) :: OrbBasis  !Gaussian basis used to expand the MOs
                                       !with centers on the entire supercell
    end type

    type, bind(c) :: op_matrix_t
        type(dyn_array_t) :: Data
        integer :: nRows
        integer :: nCols
    end type

    type, bind(c) :: SolidModel
        type(unit_cell_t) :: UnitCell
        type(lattice_t) :: Lattice
        type(supercell_t) :: Supercell
    end type

end module
