program test

REAL(8)    :: A,B
     integer :: num_args, ix
     character(len=12), dimension(:), allocatable :: args

     num_args = command_argument_count()
     allocate(args(num_args))  ! I've omitted checking the return status of the allocation 

     do ix = 1, num_args
         call get_command_argument(ix,args(ix))
         ! now parse the argument as you wish
         PRINT*, args(ix)
     end do

     PRINT*, A+B, COMMAND_ARGUMENT_COUNT()

end program test

