!
!   my_fortran_lib.f90:
!
program main

  use io_module, only:write_file, read_file
  use points_container_module
  use mesh_smoothing

  implicit none

  TYPE(points_container), pointer :: points
  character(200) :: filename
  character(15) :: str
  integer :: m, i, j

  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT() < 1)THEN

    WRITE(*,*)'ERROR, A COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
    STOP

  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,filename)
  call read_file(filename, points)
  call smooth_mesh(points, 10000, 1E-5)
  call write_file(filename, points)
  call deallocate_memory(points)

end program main
