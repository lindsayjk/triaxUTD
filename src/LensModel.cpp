#include "LensModel.h"

LensModel::LensModel(Scalar c, Scalar r200, Scalar M200, Scalar z, const Cosmology* cosmo)
	: c(c), r200(r200), M200(M200), z(z), cosmo(cosmo)
{
	if(isnan(this->r200)) LensModel::r200ToM200();
	else if(isnan(this->M200)) LensModel::M200tor200();
	else throw "Either r200 or M200 or both must be specified";
	if(!cosmo) throw "Cosmology must be specified";
}
