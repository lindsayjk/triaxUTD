#include "LensModel.h"

LensModel::LensModel()
{
}

LensModel::~LensModel()
{
}

void LensModel::setParameters(Scalar c, Scalar r200, Scalar M200, Scalar z, Scalar Dl, Scalar rhoC)
{
	this->c = c;
	this->r200 = r200;
	this->M200 = M200;
	this->z = z;
	this->Dl = Dl;
	this->rhoC = rhoC;
	rs = r200 / c;
	deltaC = (200./3.)*pow(c,3.0)/(log(1.0+c)-(c/(1.0+c)));
	rs_rhoC_deltaC = rs * rhoC * deltaC;
	if(isnan(this->r200) && !isnan(this->M200)) LensModel::M200Tor200();
	else if(isnan(this->M200) && !isnan(this->r200)) LensModel::r200toM200();
	else if(isnan(this->r200) && isnan(this->M200)) throw_line("Either r200 or M200 or both must be specified");
}
