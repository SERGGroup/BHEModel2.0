
module char_buffer

    implicit None
    public :: circular_buffer, buffer_node, allocate_buffer, deallocate_buffer, retreive_value, store_value
    private

    type buffer_node

        character(50) :: line
        TYPE (buffer_node), POINTER :: next

    end type

    type circular_buffer

        TYPE (buffer_node), POINTER :: head
        TYPE (buffer_node), POINTER :: tail

    end type

contains

    subroutine allocate_buffer(buffer, nodes)

        TYPE (circular_buffer), POINTER, intent(out) :: buffer
        TYPE (buffer_node), POINTER :: new_ptr
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

            call allocate_node(buffer)
            i = i + 1

        end do

        !connect head and tail and set head == tail
        buffer%tail%next => buffer%head
        buffer%tail => buffer%head

    end subroutine allocate_buffer

    subroutine deallocate_buffer(buffer)

        TYPE (circular_buffer), POINTER :: buffer
        TYPE (buffer_node), POINTER :: new_ptr, tail_ptr, tmp_ptr

        ! reset head position
        tail_ptr => buffer%tail
        buffer%head => tail_ptr%next

        ! disconnect tail
        nullify(new_ptr)
        tail_ptr%next => new_ptr

        ! deallocate nodes and buffer
        call deallocate_nodes(buffer%head)
        DEALLOCATE(buffer)

    end subroutine deallocate_buffer

    subroutine retreive_value(buffer, character)

        TYPE (circular_buffer), POINTER, intent(in) :: buffer
        character(50), intent(out) :: character

        character = buffer%head%line
        buffer%head => buffer%head%next

    endsubroutine retreive_value

    subroutine store_value(buffer, character)

        TYPE (circular_buffer), POINTER, intent(in) :: buffer
        character(50), intent(in) :: character

        buffer%tail%line = character
        buffer%tail => buffer%tail%next

    endsubroutine store_value

    subroutine allocate_node(buffer)

        TYPE (circular_buffer), POINTER, intent(in) :: buffer
        TYPE (buffer_node), POINTER :: new_ptr, tmp_ptr

        ! allocate new node
        ALLOCATE(new_ptr)

        ! reset tail and connect new node
        ! (new node will become the new tail, the "next" pointer
        ! of the old tail is connected to the new node)

        tmp_ptr => buffer%tail
        tmp_ptr%next => new_ptr
        buffer%tail => new_ptr

    end subroutine allocate_node

    subroutine deallocate_nodes(head_node)

        TYPE (buffer_node), POINTER :: head_node, tmp_node

        do while (ASSOCIATED(head_node))

            tmp_node => head_node%next
            DEALLOCATE(head_node)
            head_node => tmp_node

        end do


    end subroutine deallocate_nodes

end module char_buffer