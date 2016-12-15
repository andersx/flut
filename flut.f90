module flut

    implicit none

    double precision, allocatable, dimension(:) :: lut
    double precision :: lut_delta
    double precision :: lut_inv_delta
    double precision :: lut_dmin

contains

    subroutine initialize_lut(dmin, dmax, steps)

    implicit none

    double precision, intent(in) :: dmax
    double precision, intent(in) :: dmin
    integer, intent(in) :: steps

    integer :: i

    double precision :: x0

    lut_delta =  (dmax - dmin) / steps
    lut_inv_delta =  steps / (dmax - dmin)
    lut_dmin = dmin

    allocate(lut(steps))

    !$OMP PARALLEL DO
    do i = 1, steps
        lut(i) = exp(lut_dmin + (i - 1) * lut_delta)

        write (*,*) i, dmin + (i - 1) * lut_delta, lut(i)
    enddo
    !$OMP END PARALLEL DO

    end subroutine initialize_lut


    function exp_lut(x) result(e)

        implicit none

        double precision, intent(in) :: x

        double precision :: x_lut
        integer :: x_floor

        double precision :: e

        x_lut = (x - lut_dmin) * lut_inv_delta
        x_floor = floor(x_lut) + 1

        write(*,*) "lol1", x_floor, lut(x_floor)
        write(*,*) "lol2", x_floor+1, lut(x_floor+1)

        e = lut(x_floor) - (lut(x_floor + 1) - lut(x_floor)) * (x_lut - x_floor)
        print *, x_lut - x_floor

    end function exp_lut

    subroutine free_lut()

        implicit none

        deallocate(lut)
    end subroutine free_lut

end module flut


program example

    use flut, only: initialize_lut, exp_lut, free_lut

    implicit none

    double precision, parameter :: dmax = 0.0d0
    double precision, parameter :: dmin = 2.0d0
    integer, parameter :: steps = 10

    double precision :: x0, e

    call initialize_lut(dmin, dmax, steps)

    x0 = 1.115d0
    e = exp_lut(x0)
    
    write(*,*) exp(x0), e, exp(x0) / e

    call free_lut()

end program example
