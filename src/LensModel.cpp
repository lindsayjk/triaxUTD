#include "LensModel.h"

LensModel::LensModel(Scalar c, Scalar r200, Scalar M200, Scalar z, const Cosmology* cosmo)
	: c(c), r200(r200), M200(M200), z(z), cosmo(cosmo)
{
	if(!cosmo) throw "Cosmology must be specified";
	Dl = cosmo->angularDiameterDistance(z) / cosmo->h();
	rs = r200 / c;
	rhoC = cosmo->rho_c(z)*(1.E9)*cosmo->h2();
	deltaC = (200./3.)*pow(c,3.0)/(log(1.0+c)-(c/(1.0+c)));
	rs_rhoC_deltaC = rs * rhoC * deltaC;
	if(isnan(this->r200)) LensModel::r200ToM200();
	else if(isnan(this->M200)) LensModel::M200tor200();
	else throw "Either r200 or M200 or both must be specified";
}
