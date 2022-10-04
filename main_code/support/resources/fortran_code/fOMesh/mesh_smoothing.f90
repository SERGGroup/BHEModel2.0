! Created by  on 21/04/2022.

module mesh_smoothing

    use points_container_module
    implicit None

contains

    subroutine smooth_mesh(points, max_iteration, toll)

        TYPE(points_container), intent(in), pointer :: points
        integer, intent(in) :: max_iteration
        real, intent(in) :: toll

        TYPE(points_container), pointer :: new_points
        real :: residual, old_residual
        logical :: conv
        integer :: i


        call direct_memory_allocation(points%eta_dim, points%nu_dim, new_points)
        call duplicate_values(points, new_points)

        i = 0
        conv = .false.
        do while(i < max_iteration .and. .not. conv)

            call mesh_smoothing_passage(points, new_points, residual)

            residual = sqrt(residual)
            !print*, i, " -> ", residual

            if (i > 1) then

                conv = conv .or. (residual > old_residual)

            end if

            conv = conv .or. (residual < toll)
            old_residual = residual

            i = i + 1

        end do

        call deallocate_memory(new_points)

    end subroutine smooth_mesh

    subroutine mesh_smoothing_passage(points, new_points, residual)

        TYPE(points_container), intent(in), pointer :: points, new_points
        real, intent(out) :: residual

        call smooth_points(points, new_points, residual)
        call evaluate_boundary_conditions(points, new_points)
        call switch_pointers(points, new_points)

    end subroutine mesh_smoothing_passage

    subroutine smooth_points(points, new_points, residual)

        TYPE(points_container), intent(in), pointer :: points, new_points
        real, intent(out) :: residual
        real point_residual
        integer :: i, j

        i = 1
        residual = 0
        do while(i < points%eta_dim - 1)

            j = 1
            do while(j < points%nu_dim - 1)

                call smooth_point(points, new_points, i, j, point_residual)
                residual = max(residual, point_residual)
                j = j + 1

            end do
            i = i + 1

        end do

    end subroutine smooth_points

    subroutine smooth_point(points, new_points, i, j, residual)

        TYPE(points_container), intent(in), pointer :: points, new_points
        real, intent(out) :: residual
        integer, intent(in) :: i, j

        real :: dx_vert, dx_oriz, dx_tot, dy_vert, dy_oriz, dy_tot, a, b, g, smoothing
        integer :: n, j_max, n_plus, n_minus, n_right, n_left
        
        j_max = points%nu_dim

        n = i * j_max + j + 1
        n_plus = (i + 1) * j_max + j + 1
        n_minus = (i - 1) * j_max + j + 1
        n_right = i * j_max + (j + 1) + 1
        n_left = i * j_max + (j - 1) + 1

        dx_vert = points%x(n_plus) - points%x(n_minus)
        dx_oriz = points%x(n_right) - points%x(n_left)

        dy_vert = points%y(n_plus) - points%y(n_minus)
        dy_oriz = points%y(n_right) - points%y(n_left)

        dx_tot = points%x(n_plus + 1) + points%x(n_minus - 1) - points%x(n_plus - 1) - points%x(n_minus + 1)
        dy_tot = points%y(n_plus + 1) + points%y(n_minus - 1) - points%y(n_plus - 1) - points%y(n_minus + 1)

        a = 1. / 4. * (dx_oriz * dx_oriz + dy_oriz * dy_oriz)
        g = 1. / 4. * (dx_vert * dx_vert + dy_vert * dy_vert)
        b = 1. / 16. * (dx_vert * dx_oriz + dy_vert * dy_oriz)
        smoothing = 1. / (2. * (a + g + 1E-9))

        new_points%x(n) = - smoothing * (2. * b * dx_tot - a * dx_vert - g * dx_oriz)
        new_points%y(n) = - smoothing * (2. * b * dy_tot - a * dy_vert - g * dy_oriz)

        residual = (new_points%x(n) - points%x(n)) ** 2 + (new_points%y(n) - points%y(n)) ** 2

    end subroutine smooth_point

    subroutine evaluate_boundary_conditions(points, new_points)

        TYPE(points_container), intent(in), pointer :: points, new_points
        integer :: i, j, j_max, i_max, n, n_plus, n_minus, n_int
        real :: x_mean, dot_prod, delta_n

        j_max = points%nu_dim
        i_max = points%eta_dim

        ! top - bottom boundary
        ! (bottom -> fixed points, hence no modification is necessary)
        ! (top -> move the point on the line in order to make the line normal - edge of the figure are fixed)
        ! (!! to be improved freeing the tangential value of the top boundary!!)

!        i = 1
!        do while(i < i_max - 1)
!
!            n = i * j_max + j_max
!            n_plus = (i + 1) * j_max + j_max
!            n_minus = (i - 1) * j_max + j_max
!            n_int = i * j_max + j_max - 1
!
!            dot_prod_plus = (new_points%x(n_plus) - new_points%x(n))*(new_points%x(n_int) - new_points%x(n)) + (new_points%y(n_plus) - new_points%y(n))*(new_points%y(n_int) - new_points%y(n))
!            dot_prod_minus = (new_points%x(n_minus) - new_points%x(n))*(new_points%x(n_int) - new_points%x(n)) + (new_points%y(n_minus) - new_points%y(n))*(new_points%y(n_int) - new_points%y(n))
!            dot_prod_plus =
!
!            !right boundary
!            n = j
!            new_points%x(n) = x_mean
!
!            !left boundary
!            n = (i_max - 1) * j_max + j
!            new_points%x(n) = x_mean
!
!            i = i + 1
!
!        end do

        ! left - right boundary
        ! (simmetry with respect of the x axsis, y values remain unchainged)
        ! (!! to be improved as simmetry with respect to normal axis!!)

        j = 0
        do while(j < j_max)

            j = j + 1

            n_plus = 1 * j_max + j
            n_minus = (i_max - 2) * j_max + j

            x_mean = (new_points%x(n_minus) + new_points%x(n_plus)) / 2
            new_points%x(n_minus) = x_mean
            new_points%x(n_plus) = x_mean

            !right boundary
            n = j
            new_points%x(n) = x_mean

            !left boundary
            n = (i_max - 1) * j_max + j
            new_points%x(n) = x_mean

        end do

    end subroutine evaluate_boundary_conditions

end module mesh_smoothing