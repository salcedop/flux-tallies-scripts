
!===============================================================================
! FATAL_ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine fatal_error(message, error_code)

    character(*) :: message
    integer, optional :: error_code ! error code

    integer :: code      ! error code
    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message
    integer :: indent    ! length of indentation
#ifdef OPENMC_MPI
    integer :: mpi_err
#endif


    ! set default error code
    if (present(error_code)) then
      code = error_code
    else
      code = -1
    end if

    ! Write error at beginning
    write(ERROR_UNIT, fmt='(1X,A)', advance='no') 'ERROR: '

    ! Set line wrapping and indentation
    line_wrap = 80
    indent = 8

    ! Determine length of message
    length = len_trim(message)

    i_start = 0
    do
      if (length - i_start < line_wrap - indent + 1) then
        ! Remainder of message will fit on line
        write(ERROR_UNIT, fmt='(A)') message(i_start+1:length)
        exit

      else
        ! Determine last space in current line
        i_end = i_start + index(message(i_start+1:i_start+line_wrap-indent+1), &
             ' ', BACK=.true.)

        if (i_end == i_start) then
          ! This is a special case where there is no space
          i_end = i_start + line_wrap - indent + 1
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
          i_end = i_end - 1
        else
          ! Write up to last space
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
        end if

        ! Advance starting position
        i_start = i_end
        if (i_start > length) exit
      end if
    end do

#ifdef OPENMC_MPI
    ! Abort MPI
    call MPI_ABORT(mpi_intracomm, code, mpi_err)
#endif

    ! Abort program
#ifdef NO_F2008
    stop
#else
    error stop
#endif

  end subroutine fatal_error

  subroutine fatal_error_from_c(message, message_len) bind(C)
    integer(C_INT),         intent(in), value :: message_len
    character(kind=C_CHAR), intent(in)        :: message(message_len)
    character(message_len+1)                  :: message_out
    write(message_out, *) message
    call fatal_error(message_out)
  end subroutine

!===============================================================================
! WRITE_MESSAGE displays an informational message to the log file and the
! standard output stream.
!===============================================================================

  subroutine write_message(message, level)
    character(*), intent(in)           :: message ! message to write
    integer,      intent(in), optional :: level   ! verbosity level

    integer :: i_start    ! starting position
    integer :: i_end      ! ending position
    integer :: line_wrap  ! length of line
    integer :: length     ! length of message
    integer :: last_space ! index of last space (relative to start)

    ! Set length of line
    line_wrap = 80

    ! Only allow master to print to screen
    if (.not. master .and. present(level)) return

    if (.not. present(level) .or. level <= verbosity) then
      ! Determine length of message
      length = len_trim(message)

      i_start = 0
      do
        if (length - i_start < line_wrap + 1) then
          ! Remainder of message will fit on line
          write(OUTPUT_UNIT, fmt='(1X,A)') message(i_start+1:length)
          exit

        else
          ! Determine last space in current line
          last_space = index(message(i_start+1:i_start+line_wrap), &
               ' ', BACK=.true.)
          if (last_space == 0) then
            i_end = min(length + 1, i_start+line_wrap) - 1
            write(OUTPUT_UNIT, fmt='(1X,A)') message(i_start+1:i_end)
          else
            i_end = i_start + last_space
            write(OUTPUT_UNIT, fmt='(1X,A)') message(i_start+1:i_end-1)
          end if

          ! Write up to last space

        
          ! Advance starting position
          i_start = i_end
          if (i_start > length) exit
        end if
      end do
    end if

  end subroutine write_message

end module error
