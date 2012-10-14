module FortCons
    implicit none


    integer, parameter :: sp = selected_real_kind(6,37)
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: qp = selected_real_kind(33,4931)
    integer, parameter :: int32 = selected_int_kind(8)
    integer, parameter :: int64 = selected_int_kind(15)

end module FortCons
