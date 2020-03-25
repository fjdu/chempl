module ode_wrapper

USE, INTRINSIC :: ISO_C_BINDING
implicit none

contains

subroutine fdlsodes_w( &
  f, NEQ, y, t, tout, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, &
  RWORK, LRW, IWORK, LIW, jac, MF) bind(C, name="dlsodes_w")
  bind(C) f, jac
  EXTERNAL f, jac
  real(c_double), intent(inout) :: t
  real(c_double), intent(in) :: tout
  real(c_double), intent(in) :: RTOL, ATOL
  integer(c_int), intent(in) :: NEQ, ITOL, ITASK, IOPT, LRW, LIW, MF
  integer(c_int), intent(inout) :: ISTATE
  real(c_double), dimension(LRW), intent(inout) :: RWORK
  integer(c_int), dimension(LIW), intent(inout) :: IWORK
  real(c_double), dimension(NEQ), intent(inout) :: y
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
  integer(c_int), intent(in) :: mflag
  CALL XSETF(mflag)
end subroutine fxsetf_w

subroutine fxsetun_w(lun) bind(C, name="xsetun_w")
  integer(c_int), intent(in) :: lun
  CALL XSETUN(lun)
end subroutine fxsetun_w

subroutine fdsrcms_w(rsav, lrsav, isav, lisav, job) bind(C, name="dsrcms_w")
  real(c_double), dimension(lrsav), intent(inout) :: rsav
  integer(c_int), dimension(lisav), intent(inout) :: isav
  integer(c_int), intent(in) :: job, lrsav, lisav
  CALL DSRCMS(rsav, isav, job)
end subroutine fdsrcms_w

end module ode_wrapper
