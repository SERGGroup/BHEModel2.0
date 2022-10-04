! Created by  on 21/04/2022.

module real_buffer

    implicit None
    public :: real_circular_buffer, real_buffer_node, allocate_real_buffer, deallocate_real_buffer, retreive_real_value, store_real_value
    private

    type real_buffer_node

        real :: value_x, value_y
        TYPE (real_buffer_node), POINTER :: next

    end type

    type real_circular_buffer

        TYPE (real_buffer_node), POINTER :: head
        TYPE (real_buffer_node), POINTER :: tail

    end type

contains

    subroutine allocate_real_buffer(buffer, nodes)

        TYPE (real_circular_buffer), POINTER, intent(out) :: buffer
        TYPE (real_buffer_node), POINTER :: new_ptr
        integer, intent(in) :: nodes
        integer :: i

        ALLOCATE(buffer)

        !allocate first node
        !(at the beginnig it will be both head AND tail)
        ALLOCATE(new_ptr)
        buffer%tail => new_ptr
        buffer%head => new_ptr

        !allocate other nodes
        !(tail will move so that the new node will become the new tail)
        i = 1
        do while(i < nodes)

            call allocate_real_node(buffer)
            i = i + 1

        end do

        !connect head and tail and set head == tail
        buffer%tail%next => buffer%head
        buffer%tail => buffer%head

    end subroutine allocate_real_buffer

    subroutine deallocate_real_buffer(buffer)

        TYPE (real_circular_buffer), POINTER :: buffer
        TYPE (real_buffer_node), POINTER :: new_ptr, tail_ptr, tmp_ptr

        ! reset head position
        tail_ptr => buffer%tail
        buffer%head => tail_ptr%next

        ! disconnect tail
        nullify(new_ptr)
        tail_ptr%next => new_ptr

        ! deallocate nodes and buffer
        call deallocate_real_nodes(buffer%head)
        DEALLOCATE(buffer)

    end subroutine deallocate_real_buffer

    subroutine retreive_real_value(buffer, value_x, value_y)

        TYPE (real_circular_buffer), POINTER, intent(in) :: buffer
        real, intent(out) :: value_x, value_y

        value_x = buffer%head%value_x
        value_y = buffer%head%value_y
        buffer%head => buffer%head%next

    endsubroutine retreive_real_value

    subroutine store_real_value(buffer, value_x, value_y)

        TYPE (real_circular_buffer), POINTER, intent(in) :: buffer
        real, intent(in) :: value_x, value_y

        buffer%tail%value_x = value_x
        buffer%tail%value_y = value_y
        buffer%tail => buffer%tail%next

    endsubroutine store_real_value

    subroutine allocate_real_node(buffer)

        TYPE (real_circular_buffer), POINTER, intent(in) :: buffer
        TYPE (real_buffer_node), POINTER :: new_ptr, tmp_ptr

        ! allocate new node
        ALLOCATE(new_ptr)

        ! reset tail and connect new node
        ! (new node will become the new tail, the "next" pointer
        ! of the old tail is connected to the new node)

        tmp_ptr => buffer%tail
        tmp_ptr%next => new_ptr
        buffer%tail => new_ptr

    end subroutine allocate_real_node

    subroutine deallocate_real_nodes(head_node)

        TYPE (real_buffer_node), POINTER :: head_node, tmp_node

        do while (ASSOCIATED(head_node))

            tmp_node => head_node%next
            DEALLOCATE(head_node)
            head_node => tmp_node

        end do


    end subroutine deallocate_real_nodes

end module real_buffer