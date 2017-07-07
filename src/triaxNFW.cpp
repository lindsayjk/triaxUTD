#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "triaxNFW.h"
#include "mpi_logging.h"
#include "perfprof.h"

// Code for debugging integrals and logging errors and integrand values
#define CATCH_INTEGRAL_ERRORS CATCH_GSL_ERRORS

#if CATCH_INTEGRAL_ERRORS
enum IntegralType {
	IntegralNone,
	IntegralJ0,
	IntegralJ1,
	IntegralK0,
	IntegralK1,
	IntegralK2
};

struct IntegrandInfo {
	IntegralType integral_type;
	int gsl_errno;
	Scalar x2, y2, q2, r200, c, sourceSigmaC;
	Scalar a, b, theta, phi, z, Dl, rhoC;
	Scalar f, A, B, C, qX, qY;
	std::vector<Vector2> integrand_values;
};

static void reset_integrandinfo(IntegrandInfo& info, IntegralType type)
{
	info.integral_type = type;
	info.integrand_values.clear();
	info.gsl_errno = 0;
}

static IntegrandInfo LastSuccessfulJ0IntegrandInfo;
static bool HasSuccessfulJ0IntegrandInfo = false;
static IntegrandInfo LastSuccessfulJ1IntegrandInfo;
static bool HasSuccessfulJ1IntegrandInfo = false;
static IntegrandInfo LastSuccessfulK0IntegrandInfo;
static bool HasSuccessfulK0IntegrandInfo = false;
static IntegrandInfo LastSuccessfulK1IntegrandInfo;
static bool HasSuccessfulK1IntegrandInfo = false;
static IntegrandInfo LastSuccessfulK2IntegrandInfo;
static bool HasSuccessfulK2IntegrandInfo = false;

static void log_integrand_values(int gsl_errno, IntegrandInfo* info)
{
	static int log_counter = 0;

	info->gsl_errno = gsl_errno;

	// Sort the integrand values by the independent variable in ascending order
	std::sort(info->integrand_values.begin(), info->integrand_values.end(), [](const Vector2& a, const Vector2&  b) -> bool {
		return a.x<b.x;
	});

	const char* integral_name;
	switch (info->integral_type) {
	case IntegralJ0: integral_name="J0"; break;
	case IntegralJ1: integral_name="J1"; break;
	case IntegralK0: integral_name="K0"; break;
	case IntegralK1: integral_name="K1"; break;
	case IntegralK2: integral_name="K2"; break;
	default: integral_name="Unknown"; break;
	}

	char lfname[64];
	snprintf(lfname, 64, "%s.%d.integral_info", integral_name, log_counter);
	log_counter++;

	mpi_log_file lf = mpi_open_log_file(lfname);
	mpi_log(lf, "Log of %s integral", integral_name);
	mpi_log(lf, "GSL errno %d (%s)", info->gsl_errno, gsl_strerror(info->gsl_errno));
	mpi_log(lf, "x2=%g y2=%g q2=%g r200=%g c=%g sourceSigmaC=%g", info->x2, info->y2, info->q2, info->r200, info->c, info->sourceSigmaC);
	mpi_log(lf, "a=%g b=%g theta=%g phi=%g z=%g Dl=%g rhoC=%g", info->a, info->b, info->theta, info->phi, info->z, info->Dl, info->rhoC);
	mpi_log(lf, "f=%g A=%g B=%g C=%g qX=%g qY=%g", info->f, info->A, info->B, info->C, info->qX, info->qY);
	mpi_log(lf, "%d integrand values follow. Format is: u, integrand(u)", (int)info->integrand_values.size());
	for(const auto& value : info->integrand_values) {
		mpi_log(lf, "%g, %g", value.x, value.y);
	}
	mpi_close_log_file(lf);
}

static void log_last_successful_integrand_values()
{
	if(HasSuccessfulJ0IntegrandInfo) log_integrand_values(0, &LastSuccessfulJ0IntegrandInfo);
	if(HasSuccessfulJ1IntegrandInfo) log_integrand_values(0, &LastSuccessfulJ1IntegrandInfo);
	if(HasSuccessfulK0IntegrandInfo) log_integrand_values(0, &LastSuccessfulK0IntegrandInfo);
	if(HasSuccessfulK1IntegrandInfo) log_integrand_values(0, &LastSuccessfulK1IntegrandInfo);
	if(HasSuccessfulK2IntegrandInfo) log_integrand_values(0, &LastSuccessfulK2IntegrandInfo);
}

static bool integral_error_handler(int gsl_errno, void* info)
{
	if(gsl_errno==18) return false; // ignore roundoff error

	log_integrand_values(gsl_errno, static_cast<IntegrandInfo*>(info));
	log_last_successful_integrand_values();

	return true;
}

static inline void insert_integrand_info_value(const IntegrandInfo* info, Scalar u, Scalar result)
{
	Vector2 v;
	v.x = u;
	v.y = result;
	(const_cast<IntegrandInfo*>(info))->integrand_values.push_back(v);
}

static inline void save_successful_integral_info(const IntegrandInfo& info)
{
	if(info.gsl_errno!=0) return;
	if (info.integral_type == IntegralJ0) {
		if(HasSuccessfulJ0IntegrandInfo) return;
		LastSuccessfulJ0IntegrandInfo = info;
		HasSuccessfulJ0IntegrandInfo = true;
	}
	else if (info.integral_type == IntegralJ1) {
		if(HasSuccessfulJ1IntegrandInfo) return;
		LastSuccessfulJ1IntegrandInfo = info;
		HasSuccessfulJ1IntegrandInfo = true;
	}
	else if (info.integral_type == IntegralK0) {
		if(HasSuccessfulK0IntegrandInfo) return;
		LastSuccessfulK0IntegrandInfo = info;
		HasSuccessfulK0IntegrandInfo = true;
	}
	else if (info.integral_type == IntegralK1) {
		if(HasSuccessfulK1IntegrandInfo) return;
		LastSuccessfulK1IntegrandInfo = info;
		HasSuccessfulK1IntegrandInfo = true;
	}
	else if (info.integral_type == IntegralK2) {
		if(HasSuccessfulK2IntegrandInfo) return;
		LastSuccessfulK2IntegrandInfo = info;
		HasSuccessfulK2IntegrandInfo = true;
	}
}

#else
#define reset_integrandinfo(...)
#define log_integrand_values(...)
#define insert_integrand_info_value(...)
#define save_successful_integral_info(...)
#endif

#if ENABLE_PERF_PROF
DECLARE_PERF_PROF_COUNTER(J0);
DECLARE_PERF_PROF_COUNTER(J1);
DECLARE_PERF_PROF_COUNTER(K0);
DECLARE_PERF_PROF_COUNTER(K1);
DECLARE_PERF_PROF_COUNTER(K2);
#endif

// This is an alternative to GSL
static inline Scalar calc_derivative(const gsl_function* f, Scalar x, Scalar width)
{
	Scalar left = f->function(x-width/2, f->params);
	Scalar right = f->function(x+width/2, f->params);
	Scalar derivative = (right-left)/width;
	return derivative;
}

// This is what scipy was doing internally (when a and b are finite)
// a and b should be finite
static inline Scalar quad_integrate(const gsl_function* f, double a, double b, gsl_integration_workspace* w)
{
	Scalar result, abserr;
	size_t neval;
	gsl_integration_qags(f, a, b, 0, 1e-5, 500, w, &result, &abserr);
	return result;
}

// This is an alternative to GSL
// Integrate using Simpson's rule.
// Breaks the interval [from, to] into sub-intervals of length at most "width"
// and does Simpson's rule on each sub-interval.
static inline Scalar simpson_integrate(const gsl_function* f, Scalar from, Scalar to, Scalar width)
{
	Scalar signed_total_width = to-from;
	Scalar total_width = fabs(signed_total_width);
	int num_slices = ceil(total_width/width);
	Scalar result = 0.0;
	Scalar a, b, f_a, f_b, f_mid;
	b = from;
	f_b = f->function(b, f->params);
	for(int n=0;n<num_slices;n++) {
		// Integrate using Simposon's rule from a to b
		// Calculate both limits explicitly to prevent rounding errors from accumulating
		a = b;
		f_a = f_b;
		if(isnan(f_a)) {
			GSL_ERROR_VAL("simpson integration f_a is NAN", 1100, NAN);
		}
		b = from+(n+1)*signed_total_width/num_slices;
		f_b = f->function(b, f->params);
		if(isnan(f_b)) {
			GSL_ERROR_VAL("simpson integration f_b is NAN", 1101, NAN);
		}
		f_mid = f->function((a+b)/2, f->params);
		if(isnan(f_mid)) {
			GSL_ERROR_VAL("simpson integration f_mid is NAN", 1102, NAN);
		}
		result += (f_a + 4*f_mid + f_b)*(b-a)/6;
	}
	if(isnan(result)) {
		GSL_ERROR_VAL("simpson integration result is NAN", 1103, NAN);
	}
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

void triaxNFW::setParameters(Scalar c, Scalar r200, Scalar M200, Scalar a, Scalar b, Scalar theta, Scalar phi, Scalar z, Scalar Dl, Scalar rhoC)
{
	sphNFW::setParameters(c, r200, M200, z, Dl, rhoC);
	this->a = a;
	this->b = b;
	this->theta = theta;
	this->phi = phi;
	calc_fABC(theta, phi, a, b, f, A, B, C);
	inv_sqrt_f = 1/sqrt(f);
	calc_qX2_qY2(f, A, B, C, qX2, qY2);
	qX = sqrt(qX2);
	qY = sqrt(qY2);
	q2 = qY2/qX2;
	q = qY/qX;
	// Feroz and Hobson use the formula exactly as given below.
	// The Python code has atan2(B, A-C) which does not give the same angle.
	psi = atan(B/(A-C))/2;
}

triaxNFW::triaxNFW()
{
	// TODO: See if changing the workspace size has an effect on integrals
	integration_workspace = gsl_integration_workspace_alloc(1000);
}

triaxNFW::~triaxNFW()
{
	gsl_integration_workspace_free(static_cast<gsl_integration_workspace*>(integration_workspace));
}

// x2 and y2 should be prescaled by qx2
static inline Scalar calc_zeta_squared(Scalar u, Scalar x2, Scalar y2, Scalar one_minus_q2)
{
	return u*(x2+y2/(1-one_minus_q2*u));
}

// x2 and y2 should be prescaled by qx2
struct JK_integrands_params {
	triaxNFW* lens;
#if CATCH_INTEGRAL_ERRORS
	IntegrandInfo integrand_info;
#endif
	Scalar x2, y2, one_minus_q2, sourceSigmaC;
	gsl_function sph_convergence_function;
};

// The 1/sqrt(f) factor is a constant and has been pulled out of the integrands.

/*static*/ Scalar triaxNFW::J0_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar result = params->lens->calcConvergence(zeta, params->sourceSigmaC)/sqrt(1-params->one_minus_q2*u);
	insert_integrand_info_value(&params->integrand_info, u, result);
	return result;
}

/*static*/ Scalar triaxNFW::J1_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar result = params->lens->calcConvergence(zeta, params->sourceSigmaC)/pow(1-params->one_minus_q2*u, 1.5);
	insert_integrand_info_value(&params->integrand_info, u, result);
	return result;
}

/*static*/ Scalar triaxNFW::sph_convergence_function(Scalar u, const JK_integrands_params* params)
{
	return params->lens->calcConvergence(u, params->sourceSigmaC);
}

/*static*/ Scalar triaxNFW::K0_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	//begin_catch_gsl_errors("K0_integrand -> deriv sph_convergence_function");
	//gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	//end_catch_gsl_errors();
	delkappa = calc_derivative(&params->sph_convergence_function, zeta, 1e-5);
	Scalar result = u*delkappa/(zeta*sqrt(1-params->one_minus_q2*u));
	insert_integrand_info_value(&params->integrand_info, u, result);
	return result;
}

/*static*/ Scalar triaxNFW::K1_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	//begin_catch_gsl_errors("K1_integrand -> deriv sph_convergence_function");
	//gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	//end_catch_gsl_errors();
	delkappa = calc_derivative(&params->sph_convergence_function, zeta, 1e-5);
	Scalar result = u*delkappa/(zeta*pow(1-params->one_minus_q2*u, 1.5));
	insert_integrand_info_value(&params->integrand_info, u, result);
	return result;
}

/*static*/ Scalar triaxNFW::K2_integrand(Scalar u, const JK_integrands_params* params)
{
	Scalar zeta = sqrt(calc_zeta_squared(u, params->x2, params->y2, params->one_minus_q2));
	Scalar delkappa, delkappa_abserr;
	//begin_catch_gsl_errors("K2_integrand -> deriv sph_convergence_function");
	//gsl_deriv_central(&params->sph_convergence_function, zeta, 1e-5, &delkappa, &delkappa_abserr);
	//end_catch_gsl_errors();
	delkappa = calc_derivative(&params->sph_convergence_function, zeta, 1e-5);
	Scalar result = u*delkappa/(zeta*pow(1-params->one_minus_q2*u, 2.5));
	insert_integrand_info_value(&params->integrand_info, u, result);
	return result;
}

// x2 and y2 should be prescaled by qx2
inline void triaxNFW::calcJKintegrals(Scalar x2, Scalar y2, Scalar sourceSigmaC, Scalar& J0, Scalar& J1, Scalar& K0, Scalar& K1, Scalar& K2)
{
	typedef double (*integrand_fn)(double,void*);
	gsl_integration_workspace* wksp = static_cast<gsl_integration_workspace*>(integration_workspace);
	JK_integrands_params params;
	params.lens = this;
	params.x2 = x2;
	params.y2 = y2;
	params.one_minus_q2 = 1-q2;
	params.sourceSigmaC = sourceSigmaC;
	params.sph_convergence_function.function = reinterpret_cast<integrand_fn>(&triaxNFW::sph_convergence_function);
	params.sph_convergence_function.params = static_cast<void*>(&params);
#if CATCH_INTEGRAL_ERRORS
	params.integrand_info.x2 = x2;
	params.integrand_info.y2 = y2;
	params.integrand_info.q2 = q2;
	params.integrand_info.r200 = r200;
	params.integrand_info.c = c;
	params.integrand_info.sourceSigmaC = sourceSigmaC;
	params.integrand_info.f = f;
	params.integrand_info.A = A;
	params.integrand_info.B = B;
	params.integrand_info.C = C;
	params.integrand_info.qX = qX;
	params.integrand_info.qY = qY;
	params.integrand_info.a = a;
	params.integrand_info.b = b;
	params.integrand_info.theta = theta;
	params.integrand_info.phi = phi;
	params.integrand_info.z = z;
	params.integrand_info.Dl = Dl;
	params.integrand_info.rhoC = rhoC; 
#endif

	gsl_function integrand;
	integrand.params = static_cast<void*>(&params);

#if CATCH_GSL_ERRORS
	static char JKintegrals_params_str[200];
	sprintf(JKintegrals_params_str, "calcJKintegrals x2=%f y2=%f q2=%f r200=%f c=%f sourceSigmaC=%f", x2, y2, q2, r200, c, sourceSigmaC);
	begin_catch_gsl_errors(JKintegrals_params_str);
#endif

	START_PERF_PROF(J0);
	reset_integrandinfo(params.integrand_info, IntegralJ0);
	integrand.function = reinterpret_cast<integrand_fn>(&triaxNFW::J0_integrand);
	begin_catch_gsl_errors("J0", 0, integral_error_handler, &params.integrand_info);
	J0 = quad_integrate(&integrand, 0, 1, wksp) * inv_sqrt_f;
	end_catch_gsl_errors();
	save_successful_integral_info(params.integrand_info);
	END_PERF_PROF(J0);
	ACCUM_PERF_PROF_DURATION(J0, J0);

	START_PERF_PROF(J1);
	reset_integrandinfo(params.integrand_info, IntegralJ1);
	integrand.function = reinterpret_cast<integrand_fn>(&triaxNFW::J1_integrand);
	begin_catch_gsl_errors("J1", 0, integral_error_handler, &params.integrand_info);
	J1 = quad_integrate(&integrand, 0, 1, wksp) * inv_sqrt_f;
	end_catch_gsl_errors();
	save_successful_integral_info(params.integrand_info);
	END_PERF_PROF(J1);
	ACCUM_PERF_PROF_DURATION(J1, J1);

	START_PERF_PROF(K0);
	reset_integrandinfo(params.integrand_info, IntegralK0);
	integrand.function = reinterpret_cast<integrand_fn>(&triaxNFW::K0_integrand);
	begin_catch_gsl_errors("K0", 0, integral_error_handler, &params.integrand_info);
//	K0 = quad_integrate(&integrand, 0, 1, wksp) * inv_sqrt_f;
	K0 = simpson_integrate(&integrand, 1e-3, 1, 1./10) * inv_sqrt_f;
	end_catch_gsl_errors();
	save_successful_integral_info(params.integrand_info);
	END_PERF_PROF(K0);
	ACCUM_PERF_PROF_DURATION(K0, K0);

	START_PERF_PROF(K1);
	reset_integrandinfo(params.integrand_info, IntegralK1);
	integrand.function = reinterpret_cast<integrand_fn>(&triaxNFW::K1_integrand);
	begin_catch_gsl_errors("K1", 0, integral_error_handler, &params.integrand_info);
//	K1 = quad_integrate(&integrand, 0, 1, wksp) * inv_sqrt_f;
	K1 = simpson_integrate(&integrand, 1e-3, 1, 1./10) * inv_sqrt_f;
	end_catch_gsl_errors();
	save_successful_integral_info(params.integrand_info);
	END_PERF_PROF(K1);
	ACCUM_PERF_PROF_DURATION(K1, K1);

	START_PERF_PROF(K2);
	reset_integrandinfo(params.integrand_info, IntegralK2);
	integrand.function = reinterpret_cast<integrand_fn>(&triaxNFW::K2_integrand);
	begin_catch_gsl_errors("K2", 0, integral_error_handler, &params.integrand_info);
//	K2 = quad_integrate(&integrand, 0, 1, wksp) * inv_sqrt_f;
	K2 = simpson_integrate(&integrand, 1e-3, 1, 1./10) * inv_sqrt_f;
	end_catch_gsl_errors();
	save_successful_integral_info(params.integrand_info);
	END_PERF_PROF(K2);
	ACCUM_PERF_PROF_DURATION(K2, K2);

	end_catch_gsl_errors();
}

// x and y should be prescaled by qx
inline void triaxNFW::calcPotential2ndDerivatives(Scalar x, Scalar y, Scalar sourceSigmaC, Scalar& pot_xx, Scalar& pot_xy, Scalar& pot_yy)
{
	Scalar x2 = x*x;
	Scalar y2 = y*y;
	Scalar J0, J1, K0, K1, K2;
	calcJKintegrals(x2, y2, sourceSigmaC, J0, J1, K0, K1, K2);
	pot_xx = q*(x2*K0 + J0);
	pot_yy = q*(y2*K2 + J1);
	pot_xy = q*x*y*K1;
}

// Calculated convergence and shear in the intermediate coordinate system
// x and y should be prescaled by qx
inline void triaxNFW::calcIntermediateConvergenceShear(Scalar x, Scalar y, Scalar sourceSigmaC, Scalar& kappa, Scalar& gamma1_intermediate, Scalar& gamma2_intermediate)
{
	Scalar pot_xx, pot_xy, pot_yy;
	calcPotential2ndDerivatives(x, y, sourceSigmaC, pot_xx, pot_xy, pot_yy);
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
void triaxNFW::calcConvergenceShear(Vector2Array1DRef coord_list, ScalarArray1DRef sourceSigmaC_list, ScalarArray1DRef kappa_out, ScalarArray1DRef gamma1_out, ScalarArray1DRef gamma2_out)
{
	START_PERF_PROF(calcConvergenceShear);
	DECLARE_PERF_PROF_COUNTER(calcIntermediateConvergenceShear);
	RESET_PERF_PROF_COUNTER(J0);
	RESET_PERF_PROF_COUNTER(J1);
	RESET_PERF_PROF_COUNTER(K0);
	RESET_PERF_PROF_COUNTER(K1);
	RESET_PERF_PROF_COUNTER(K2);
	int num_coords = coord_list->getLen();
	Vector2* v = coord_list->v;
	Scalar* p_sourceSigmaC = sourceSigmaC_list->v;
	Scalar* p_kappa_out = kappa_out->v;
	Scalar* p_gamma1_out = gamma1_out->v;
	Scalar* p_gamma2_out = gamma2_out->v;

	Scalar sin_2psi = sin(2*psi);
	Scalar cos_2psi = cos(2*psi);

	int n;
	for (n = 0; n < num_coords; n++, v++, p_sourceSigmaC++, p_kappa_out++, p_gamma1_out++, p_gamma2_out++) {
		Scalar x = v->x/qX;
		Scalar y = v->y/qX;

		// transform x and y to intermediate coordinate system by rotating through -2*psi
		Scalar xi = x*cos_2psi + y*sin_2psi;
		Scalar yi = -x*sin_2psi + y*cos_2psi;

		Scalar kappa, gamma1_intermediate, gamma2_intermediate;
		START_PERF_PROF(calcIntermediateConvergenceShear);
		calcIntermediateConvergenceShear(xi, yi, *p_sourceSigmaC, kappa, gamma1_intermediate, gamma2_intermediate);
		END_PERF_PROF(calcIntermediateConvergenceShear);
		ACCUM_PERF_PROF_DURATION(calcIntermediateConvergenceShear, calcIntermediateConvergenceShear);

		Scalar gamma1, gamma2;
		// transform gamma out of the intermediate coordinate system by rotating through 2*psi
		gamma1 = gamma1_intermediate*cos_2psi - gamma2_intermediate*sin_2psi;
		gamma2 = gamma1_intermediate*sin_2psi + gamma2_intermediate*cos_2psi;
		
		*p_kappa_out=kappa;
		*p_gamma1_out=gamma1;
		*p_gamma2_out=gamma2;
	}
	END_PERF_PROF(calcConvergenceShear);
	if(ENABLE_PERF_PROF) {
		mpi_log(NULL, "1 iteration of triaxNFW::calcConvergenceShear took %lld ns, %d iterations of calcIntermediateConvergenceShear took %lld ns", GET_PERF_PROF_DURATION(calcConvergenceShear), n, PERF_PROF_COUNTER_NS(calcIntermediateConvergenceShear));
		mpi_log(NULL, "\tAll iterations of J0,J1,K0,K1,K2 integrals took %lld %lld %lld %lld %lld (ns)",
			PERF_PROF_COUNTER_NS(J0),
			PERF_PROF_COUNTER_NS(J1),
			PERF_PROF_COUNTER_NS(K0),
			PERF_PROF_COUNTER_NS(K1),
			PERF_PROF_COUNTER_NS(K2)
		);
	}
}
