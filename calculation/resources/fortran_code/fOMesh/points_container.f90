! Created by  on 21/04/2022.

module points_container_module

    implicit None
    public :: points_container, filename_memory_allocation, direct_memory_allocation, deallocate_memory, switch_pointers, duplicate_values
    private

    type points_container

        real, dimension(:), pointer :: x, y
        integer :: eta_dim, nu_dim, n_tot

    end type

contains

    subroutine filename_memory_allocation(filename, points)

        TYPE(points_container), intent(out), pointer :: points
        character(len=*), intent(in) :: filename

        character(30) :: other_str, final_str
        real, dimension(:), pointer :: x, y
        integer :: eta_dim, nu_dim

        call split_string(filename, '-', other_str, eta_dim)
        call split_string(other_str, '-', final_str, nu_dim)

        allocate(points)
        allocate(x(eta_dim * nu_dim))
        allocate(y(eta_dim * nu_dim))

        points%x => x
        points%y => y

        points%eta_dim = eta_dim
        points%nu_dim = nu_dim
        points%n_tot = eta_dim * nu_dim

    end subroutine filename_memory_allocation

    subroutine direct_memory_allocation(eta_dim, nu_dim, points)

        TYPE(points_container), intent(out), pointer :: points
        integer, intent(in)  :: eta_dim, nu_dim
        real, dimension(:), pointer :: x, y

        allocate(points)
        allocate(x(eta_dim * nu_dim))
        allocate(y(eta_dim * nu_dim))

        points%x => x
        points%y => y

        points%eta_dim = eta_dim
        points%nu_dim = nu_dim
        points%n_tot = eta_dim * nu_dim

    end subroutine direct_memory_allocation

    subroutine deallocate_memory(points)

        TYPE(points_container), pointer :: points

        deallocate(points%x)
        deallocate(points%y)
        deallocate(points)

    end subroutine deallocate_memory

    subroutine split_string(instring, delim, other_str, value)

        character(len=*), intent(in) :: instring
        character(30), intent(out) :: other_str
        integer, intent(out) :: value

        character(1)  :: delim
        integer :: index

        index = scan(instring,delim)
        other_str = instring(index + 1:)
        read(instring(1:index-1), *) value

    end subroutine split_string

    subroutine switch_pointers(points, new_points)

        ! switch the pointers allow to update the calculations

        TYPE(points_container), intent(in), pointer :: points, new_points
        real, dimension(:), pointer :: tmp_x, tmp_y

        tmp_x => points%x
        tmp_y => points%y

        points%x => new_points%x
        points%y => new_points%y

        new_points%x => tmp_x
        new_points%y => tmp_y

    end subroutine switch_pointers

    subroutine duplicate_values(points, new_points)

        ! copy the vaalues written in points into new_points

        TYPE(points_container), intent(in), pointer :: points, new_points
        real, dimension(:), pointer :: tmp_x, tmp_y
        integer :: n

        n = 0
        do while(n < points%n_tot)

            n = n + 1
            new_points%x(n) = points%x(n)
            new_points%y(n) = points%y(n)

        end do

    end subroutine duplicate_values


end module points_container_module