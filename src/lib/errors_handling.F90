module error_mod
    implicit none
    contains

    subroutine stop_all(sub_name,error_msg)
        != Stop calculation due to an error.
        != Exit with code 999.
        !=
        != In:
        !=    sub_name:  calling subroutine name.
        !=    error_msg: error message.

        implicit none
        character(*), intent(in) :: sub_name,error_msg

        ! It seems that giving STOP a string is far more portable.
        ! MPI_Abort requires an integer though.
        character(3), parameter :: error_str='999'

        write (6,'(/a7)') 'ERROR.'
        write (6,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
        write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
        write (6,'(a11)') 'EXITING...'
        stop '999'
        return
    end subroutine stop_all

end module error_mod

