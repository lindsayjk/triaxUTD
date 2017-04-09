#include <gsl/gsl_integration.h>
#include "triaxNFW.h"

// This is what scipy was doing internally (when a and b are finite)
// a and b should be finite
static inline Scalar quad_integrate(const gsl_function* f, double a, double b, gsl_integration_workspace* w)
{
	Scalar result, abserr;
	gsl_integration_qags(f, a, b, 1.49e-8, 1.49e-8, 50, w, &result, &abserr);
	return result;
}

static inline void calc_fABC(Scalar theta, Scalar phi, Scalar a, Scalar b, Scalar& f, Scalar& A, Scalar& B, Scalar& C) // c = 1
{
	Scalar sin_theta = sin(theta);
	Scalar cos_theta = cos(theta);
	Scalar sin_phi = sin(phi);
	Scalar cos_phi = cos(phi);
	Scalar sin2_theta = sin_theta*sin_theta;
	Scalar cos2_theta = cos_theta*cos_theta;
	Scalar sin2_phi = sin_phi*sin_phi;
	Scalar cos2_phi = cos_phi*cos_phi;
	Scalar sin_2phi = sin(2*phi);
	Scalar inv_a2 = 1.0/(a*a);
	Scalar inv_b2 = 1.0/(b*b);
	f = sin2_theta*(cos2_phi*inv_a2 + sin2_phi*inv_b2) + cos2_theta;
	A = cos2_theta*(sin2_phi*inv_a2 + cos2_phi*inv_b2) + sin2_theta*inv_a2*inv_b2;
	B = cos_theta*sin_2phi*(inv_a2-inv_b2);
	C = sin2_phi*inv_b2 + cos2_phi*inv_a2;
}

static inline void calc_qX2_qY2(Scalar f, Scalar A, Scalar B, Scalar C, Scalar& qX2, Scalar& qY2)
{
	Scalar a = A + C;
	Scalar b = sqrt((A-C)*(A-C)+B*B);
	qX2 = 2*f/(a-b);
	qY2 = 2*f/(a+b);
}

triaxNFW::triaxNFW(Scalar c, Scalar r200, Scalar M200, Scalar a, Scalar b, Scalar theta, Scalar phi, Scalar z, const Cosmology* cosmo)
	: sphNFW(c, r200, M200, z, cosmo), a(a), b(b), theta(theta), phi(phi)
{
	// TODO: See if changing the workspace size has an effect on integrals
	integration_workspace = gsl_integration_workspace_alloc(1000);
	calc_fABC(theta, phi, a, b, f, A, B, C);
	inv_sqrt_f = 1/sqrt(f);
	calc_qX2_qY2(f, A, B, C, qX2, qY2);
	qX = sqrt(qX2);
	qY = sqrt(qY2);
	q2 = qY2/qX2;
	q = qY/qX;
	// DISCUSSION:
	// Feroz and Hobson use the formula exactly as given below.
	// The Python code has atan2(B, A-C) which does not give the same angle.
	// Is the Python code deliberately based on a different definition of psi or was it a mistake?
	psi = atan(B/(A-C))/2;
}

triaxNFW::~triaxNFW()
{
	gsl_integration_workspace_free(integration_workspace);
}

// x2 and y2 should be prescaled by qx2
static inline Scalar calc_zeta_squared(Scalar u, Scalar x2, Scalar y2, Scalar one_minus_q2)
{
	return u*(x2+y2/(1-one_minus_q2*u));
}

// x2 and y2 should be prescaled by qx2
struct JK_integrands_params {
	triaxNFW* lens;
	Scalar x2, y2, one_minus_q2, zs;
	gsl_function sph_convergence_function;
};

// The 1/sqrt(f) factor is a constant and has been pulled out of the integrands.

static Scalar triaxNFW::J0_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	return params->lens->calcConvergence(zeta, params->zs)/sqrt(1-params->one_minus_q2*u);
}

static Scalar triaxNFW::J1_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	return params->lens->calcConvergence(zeta, params->zs)/pow(1-params->one_minus_q2*u, 1.5);
}

static Scalar triaxNFW::sph_convergence_function(Scalar u, const JK_integrands_params* params)
{
	return params->lens->calcConvergence(u, params->zs);
}

static Scalar triaxNFW::K0_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	return u*delkappa/sqrt(zeta*(1-params->one_minus_q2)*u);
}

static Scalar triaxNFW::K1_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	return u*delkappa/pow(zeta*(1-params->one_minus_q2)*u, 1.5);
}

static Scalar triaxNFW::K2_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	return u*delkappa/pow(zeta*(1-params->one_minus_q2)*u, 2.5);
}

// x2 and y2 should be prescaled by qx2
inline void triaxNFW::calcJKintegrals(Scalar x2, Scalar y2, Scalar zs, Scalar& J0, Scalar& J1, Scalar& K0, Scalar& K1, Scalar& K2)
{
	JK_integrands_params params;
	params.lens = this;
	params.x2 = x2;
	params.y2 = y2;
	params.one_minus_q2 = 1-q2;
	params.zs = zs;
	params.sph_convergence_function.function = reinterpret_cast<gsl_function>(&triaxNFW::sph_convergence_function);
	params.sph_convergence_function.params = static_cast<void*>(&params);

	gsl_function integrand;
	integrand.params = static_cast<void*>(&params);

	integrand.function = reinterpret_cast<gsl_function>(&triaxNFW::J0_integrand);
	J0 = quad_integrate(&integrand, 0, 1, integration_workspace) * inv_sqrt_f;

	integrand.function = reinterpret_cast<gsl_function>(&triaxNFW::J1_integrand);
	J1 = quad_integrate(&integrand, 0, 1, integration_workspace) * inv_sqrt_f;

	integrand.function = reinterpret_cast<gsl_function>(&triaxNFW::K0_integrand);
	K0 = quad_integrate(&integrand, 0, 1, integration_workspace) * inv_sqrt_f;

	integrand.function = reinterpret_cast<gsl_function>(&triaxNFW::K1_integrand);
	K1 = quad_integrate(&integrand, 0, 1, integration_workspace) * inv_sqrt_f;

	integrand.function = reinterpret_cast<gsl_function>(&triaxNFW::K2_integrand);
	K2 = quad_integrate(&integrand, 0, 1, integration_workspace) * inv_sqrt_f;
}

// x and y should be prescaled by qx
inline void triaxNFW::calcPotential2ndDerivatives(Scalar x, Scalar y, Scalar zs, Scalar& pot_xx, Scalar& pot_xy, Scalar& pot_yy)
{
	Scalar x2 = x*x;
	Scalar y2 = y*y;
	Scalar J0, J1, K0, K1, K2;
	calcJKintegrals(x2, y2, zs, J0, J1, K0, K1, K2);
	pot_xx = q*(x2*K0 + J0);
	pot_yy = q*(y2*K2 + J1);
	pot_xy = q*x*y*K1;
}

// Calculated convergence and shear in the intermediate coordinate system
// x and y should be prescaled by qx
inline void triaxNFW::calcIntermediateConvergenceShear(Scalar x, Scalar y, Scalar zs, Scalar& kappa, Scalar& gamma1_intermediate, Scalar& gamma2_intermediate)
{
	Scalar pot_xx, pot_xy, pot_yy;
	calcPotential2ndDerivatives(x, y, zs, pot_xx, pot_xy, pot_yy);
	kappa = (pot_xx+pot_yy)/2;
	gamma1_intermediate = (pot_xx-pot_yy)/2;
	gamma2_intermediate = pot_xy;
}

// DISCUSSION:
// The original Python code does not seem to do a translation or rotation from x,y to the intermediate coordinate system, only scaling.
// Feroz and Hobson mention translation and rotation but only provide the equation for rotating the resulting shear back to the original coordinate system.
// I have assumed the cluster (lens) location is always at x=0, y=0 w.r.t. to the observer so the translation is not needed. Is this a valid assumption?
// However, I have added a rotation from the scaled x,y to the intermediate coordinate system by rotation through -2*psi since it is the opposite of the shear being transformed by rotating through 2*psi as given in the paper.
// Instead of calculating the intermediate shear angle/magnitude and adding 2*psi to it, I did the rotation directly by a rotation matrix. This shouldn't change the result but means we don't need to calculate arctan.
void triaxNFW::calcConvergenceShear(Vector2Array1DRef coord_list, Scalar z_source, ScalarArray1DRef kappa_out, ScalarArray1DRef gamma1_out, ScalarArray1DRef gamma2_out)
{
	int num_coords = coord_list->getLen();
	Vector2* v = coord_list->v;
	Scalar* p_kappa_out = kappa_out->v;
	Scalar* p_gamma1_out = gamma1_out->v;
	Scalar* p_gamma2_out = gamma2_out->v;

	Scalar sin_2psi = sin(2*psi);
	Scalar cos_2psi = cos(2*psi);

	for (int n = 0; n < num_coords; n++, v++, p_kappa_out++, p_gamma1_out++, p_gamma2_out++) {
		Scalar x = v->x/qX;
		Scalar y = v->y/qX;

		// transform x and y to intermediate coordinate system by rotating through -2*psi
		Scalar xi = x*cos_2psi + y*sin_2psi;
		Scalar yi = -x*sin_2psi + y*cos_2psi;

		Scalar kappa, gamma1_intermediate, gamma2_intermediate;
		calcIntermediateConvergenceShear(xi, yi, z_source, kappa, gamma1_intermediate, gamma2_intermediate);

		Scalar gamma1, gamma2;
		// transform gamma out of the intermediate coordinate system by rotating through 2*psi
		gamma1 = gamma1_intermediate*cos_2psi - gamma2_intermediate*sin_2psi;
		gamma2 = gamma1_intermediate*sin_2psi + gamma2_intermediate*cos_2psi;
		
		*p_kappa_out=kappa;
		*p_gamma1_out=gamma1;
		*p_gamma2_out=gamma2;
	}
}
