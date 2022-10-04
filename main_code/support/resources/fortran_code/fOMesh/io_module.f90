
module io_module

    use char_buffer
    use real_buffer
    use points_container_module

    implicit None
    public :: write_file, read_file
    private

contains

    subroutine write_file(filename, points)

        TYPE(points_container), intent(in), pointer :: points
        character(len=*), intent(in) :: filename
        integer :: i, iu

        open(newunit=iu, file=filename, action='write', status='replace')

        do i=1,points%n_tot

            write(iu, "(E16.8,A,E16.8)"), points%x(i),';', points%y(i)

        enddo

        close(iu)

    end subroutine write_file

    subroutine read_file(filename, points)

        TYPE(points_container), intent(out), pointer :: points
        character(len=*), intent(in) :: filename

        integer :: iu, ierr, buffersize, maxbuffersize, delta_n_max
        TYPE (circular_buffer), POINTER :: buffer
        TYPE (real_circular_buffer), POINTER :: values_buffer
        integer, pointer :: n

        maxbuffersize = 100

        call filename_memory_allocation(filename, points)
        call allocate_buffer(buffer, maxbuffersize + 50)
        call allocate_real_buffer(values_buffer, maxbuffersize + 50)
        allocate(n)

        n = 0
        ierr = 0
        open(newunit=iu, file=filename, action='read', status='old')

        do while(n < points%n_tot .and. ierr == 0)

            delta_n_max = points%n_tot - n
            call store_lines_in_buffer(buffer, iu, delta_n_max, maxbuffersize, buffersize, ierr)
            call store_values_in_buffer(buffer, values_buffer, buffersize)
            call write_values_in_pointer(values_buffer, buffersize, n, points)

        enddo

        call deallocate_real_buffer(values_buffer)
        call deallocate_buffer(buffer)
        deallocate(n)
        close(iu)

    end subroutine read_file

    subroutine store_lines_in_buffer(buffer, iu, delta_n_max, maxbuffersize, buffersize, ierr)

        TYPE (circular_buffer), POINTER, intent(in) :: buffer
        integer, intent(in) :: iu, maxbuffersize, delta_n_max
        integer, intent(out) :: ierr, buffersize
        character(50) :: new_line
        integer :: i

        i = 0
        buffersize = min(maxbuffersize, delta_n_max)

        do while(i < buffersize .and. ierr == 0)

            read(iu, '(A)', iostat=ierr), new_line
            wait(iu)
            call store_value(buffer, new_line)
            i = i + 1

        end do

        buffersize = i

    end subroutine store_lines_in_buffer

    subroutine store_values_in_buffer(buffer, values_buffer, buffersize)

        TYPE (circular_buffer), POINTER, intent(in) :: buffer
        TYPE (real_circular_buffer), POINTER, intent(in) :: values_buffer
        integer, intent(in) :: buffersize
        character(50) :: new_line
        real :: x_new, y_new
        character(3) :: tmp
        integer :: i

        i = 0
        do while(i < buffersize)

            i = i + 1
            call retreive_value(buffer, new_line)
            read(new_line, "(E16.8,A,E16.8)"), x_new, tmp, y_new
            call store_real_value(values_buffer, x_new, y_new)

        enddo

    end subroutine store_values_in_buffer

    subroutine write_values_in_pointer(buffer, buffersize, n, points)

        TYPE (real_circular_buffer), POINTER, intent(in) :: buffer
        TYPE(points_container), intent(in), pointer :: points
        integer, intent(out), pointer :: n
        integer, intent(in) :: buffersize
        real :: x_new, y_new
        integer :: i

        i = 0
        do while(i < buffersize)

            n = n + 1
            i = i + 1
            call retreive_real_value(buffer, x_new, y_new)
            points%x(n) = x_new
            points%y(n) = y_new

        enddo

    end subroutine write_values_in_pointer

end module io_module