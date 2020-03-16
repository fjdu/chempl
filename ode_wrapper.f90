module ode_wrapper

USE, INTRINSIC :: ISO_C_BINDING
implicit none

contains

subroutine fdlsodes_w( &
  f, NEQ, y, t, tout, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, &
  RWORK, LRW, IWORK, LIW, jac, MF) bind(C, name="dlsodes_w")
  bind(C) f, jac
  EXTERNAL f, jac
  double precision, intent(inout) :: t
  double precision, intent(in) :: tout
  double precision, intent(in) :: RTOL, ATOL
  integer, intent(in) :: NEQ, ITOL, ITASK, IOPT, LRW, LIW, MF
  integer, intent(inout) :: ISTATE
  double precision, dimension(LRW), intent(inout) :: RWORK
  integer, dimension(LIW), intent(inout) :: IWORK
  double precision, dimension(NEQ), intent(inout) :: y
  !integer i
  !write(*, "(A)") "Going to call DLSODES."
  !write(*,*) NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
  !do i=1, NEQ
  !  write(*,*) i, y(i)
  !end do
  call DLSODES (f, NEQ, y, t, tout, ITOL, RTOL, ATOL, &
        ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jac, MF)
  !write(*, "(A)") "Finish calling DLSODES."
end subroutine fdlsodes_w


subroutine fxsetf_w(mflag) bind(C, name="xsetf_w")
  integer, intent(in) :: mflag
  CALL XSETF(mflag)
end subroutine fxsetf_w

subroutine fxsetun_w(lun) bind(C, name="xsetun_w")
  integer, intent(in) :: lun
  CALL XSETUN(lun)
end subroutine fxsetun_w

end module ode_wrapper
