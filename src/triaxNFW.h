#pragma once
#include "sphNFW.h"

class triaxNFW : public sphNFW {
public:
	// Either r200 or M200 may be NAN but not both.
	triaxNFW(Scalar c, Scalar r200, Scalar M200, Scalar a, Scalar b, Scalar theta, Scalar phi, Scalar z, const Cosmology* cosmo);
	~triaxNFW();

	void calcConvergenceShear(Vector2Array1DRef coord_list, Scalar sourceSigmaC, ScalarArray1DRef kappa_out, ScalarArray1DRef gamma1_out, ScalarArray1DRef gamma2_out);

protected:
	Scalar a, b, theta, phi;
	Scalar f, A, B, C;
	Scalar inv_sqrt_f; // 1/sqrt(f)
	Scalar q, qX, qY;
	Scalar q2, qX2, qY2; // q, qX, and qY squared
	Scalar psi;

	struct gsl_integration_workspace* integration_workspace;

	static Scalar J0_integrand(Scalar u, const struct JK_integrands_params*);
	static Scalar J1_integrand(Scalar u, const struct JK_integrands_params*);
	static Scalar K0_integrand(Scalar u, const struct JK_integrands_params*);
	static Scalar K1_integrand(Scalar u, const struct JK_integrands_params*);
	static Scalar K2_integrand(Scalar u, const struct JK_integrands_params*);
	static Scalar sph_convergence_function(Scalar u, const struct JK_integrands_params*);

	inline void calcJKintegrals(Scalar x2, Scalar y2, Scalar sourceSigmaC, Scalar& J0, Scalar& J1, Scalar& K0, Scalar& K1, Scalar& K2);
	inline void calcPotential2ndDerivatives(Scalar x, Scalar y, Scalar sourceSigmaC, Scalar& pot_xx, Scalar& pot_xy, Scalar& pot_yy);
	inline void calcIntermediateConvergenceShear(Scalar x, Scalar y, Scalar sourceSigmaC, Scalar& kappa_intermediate, Scalar& gamma1_intermediate, Scalar& gamma2_intermediate);
};
