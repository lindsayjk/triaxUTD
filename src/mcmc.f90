module triaxUTD

use iso_c_binding
implicit none

interface

  subroutine setup() bind(c, name="triaxUTD_setup")
  end subroutine

  function lnlikelihood(c, r200, a, b, phi, theta) bind(c, name="triaxUTD_lnlikelihood")
	real(c_double) :: lnlikelihood
	real(c_double), value :: c
	real(c_double), value :: r200
	real(c_double), value :: a
	real(c_double), value :: b
	real(c_double), value :: phi
	real(c_double), value :: theta
  end function

end interface

end module triaxUTD
