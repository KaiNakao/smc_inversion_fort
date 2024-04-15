module type_mat
    implicit none
    type mat
        double precision, allocatable :: body(:,:)
        integer :: nrow, ncol
    end type
end module type_mat