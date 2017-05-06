#include "sphNFW.h"

sphNFW::sphNFW()
{
}

sphNFW::~sphNFW()
{
}

static inline Scalar calcSurfaceProfile_gt_1(Scalar x)
{
//	f = (1/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5));
	Scalar x2m1 = x*x - 1;
	return (1/x2m1)*(1 - 2*atan(sqrt((x - 1)/(x + 1)))/sqrt(x2m1));
}

static inline Scalar calcSurfaceProfile_lt_1(Scalar x)
{
//	f = (1/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5)));
	Scalar x2m1 = x*x - 1;
	return (1/x2m1)*(1 - 2*atanh(sqrt((1 - x)/(x + 1)))/sqrt(-x2m1));
}

static inline Scalar calcShear_gt_1(Scalar x)
{
//	f = (8*np.arctan(((x-1)/(1+x))**0.5)/(x**2*(x**2-1)**0.5)) + (4*np.log(x/2.)/x**2) - (2./(x**2-1)) + (4*np.arctan(((x-1)/(1+x))**0.5)/((x**2-1)**1.5));
	Scalar x2 = x*x;
	Scalar x2m1 = x2-1;
	Scalar xm1 = x-1;
	Scalar xp1 = 1+x;
	Scalar atanr = atan(sqrt(xm1/xp1));
	return (8*atanr/(x2*sqrt(x2m1))) + (4*log(x/2.0)/x2) - (2.0/x2m1) + (4*atanr/pow(x2m1, 1.5));

}

static inline Scalar calcShear_lt_1(Scalar x)
{
//	f = (8*np.arctanh(((1-x)/(1+x))**0.5)/(x**2*(1-x**2)**0.5)) + (4*np.log(x/2.)/x**2) - (2./(x**2-1)) + (4*np.arctanh(((1-x)/(1+x))**0.5)/((x**2-1)*(1-x**2)**0.5));
	Scalar x2 = x*x;
	Scalar x2m1 = x2-1;
	Scalar _1mx = 1-x;
	Scalar xp1 = 1+x;
	Scalar atanhr = atanh(sqrt(_1mx/xp1));
	Scalar sqrt_1mx2 = sqrt(-x2m1);
	return (8*atanhr/(x2*sqrt_1mx2)) + (4*log(x/2.0)/x2) - (2.0/x2m1) + (4*atanhr/(x2m1*sqrt_1mx2));
}

Scalar sphNFW::calcScaledSurfaceProfile(Scalar r, Scalar scale) const
{
	Scalar rs = get_rs();
	Scalar x = r/rs;
	if(x<1) return scale*calcSurfaceProfile_lt_1(x);
	else if(x>1) return scale*calcSurfaceProfile_gt_1(x);
	else return scale/3;
}

void sphNFW::calcScaledSurfaceProfile(ScalarArray1DRef r_array, Scalar scale, ScalarArray1DRef out_array) const
{
	Scalar* r = r_array->v;
	Scalar* out = out_array->v;
	Scalar rs = get_rs();
	int len = out_array->getLen();
	for (int n = 0; n < len; n++, out++, r++) {
		Scalar x = (*r)/rs;
		if(x<1) *out = scale*calcSurfaceProfile_lt_1(x);
		else if(x>1) *out = scale*calcSurfaceProfile_gt_1(x);
		else *out = scale/3;
	}
}

void sphNFW::calcScaledSurfaceProfile(ScalarArray2DRef r_array, Scalar scale, ScalarArray2DRef out_array) const
{
	Scalar* r = r_array->v;
	Scalar* out = out_array->v;
	Scalar rs = get_rs();
	int len = out_array->getWidth()*out_array->getHeight();
	for (int n = 0; n < len; n++, out++, r++) {
		Scalar x = (*r)/rs;
		if(x<1) *out = scale*calcSurfaceProfile_lt_1(x);
		else if(x>1) *out = scale*calcSurfaceProfile_gt_1(x);
		else *out = scale/3;
	}
}

Scalar sphNFW::calcScaledShear(Scalar r, Scalar scale) const
{
	Scalar rs = get_rs();
	Scalar x = r/rs;
	if(x<1) return scale*calcShear_lt_1(x);
	else if(x>1) return scale*calcShear_gt_1(x);
	else return scale*0.5607446110935520956644048475006270610313327958923123; // scale*(10/3+4*log(0.5))
}

void sphNFW::calcScaledShear(ScalarArray1DRef r_array, Scalar scale, ScalarArray1DRef out_array) const
{
	Scalar* r = r_array->v;
	Scalar* out = out_array->v;
	Scalar rs = get_rs();
	int len = out_array->getLen();
	for (int n = 0; n < len; n++, out++, r++) {
		Scalar x = (*r)/rs;
		if(x<1) *out = scale*calcShear_lt_1(x);
		else if(x>1) *out = scale*calcShear_gt_1(x);
		else *out = scale*0.5607446110935520956644048475006270610313327958923123; // scale*(10/3+4*log(0.5))
	}
}

void sphNFW::calcScaledShear(ScalarArray2DRef r_array, Scalar scale, ScalarArray2DRef out_array) const
{
	Scalar* r = r_array->v;
	Scalar* out = out_array->v;
	Scalar rs = get_rs();
	int len = out_array->getWidth()*out_array->getHeight();
	for (int n = 0; n < len; n++, out++, r++) {
		Scalar x = (*r)/rs;
		if(x<1) *out = scale*calcShear_lt_1(x);
		else if(x>1) *out = scale*calcShear_gt_1(x);
		else *out = scale*0.5607446110935520956644048475006270610313327958923123; // scale*(10/3+4*log(0.5))
	}
}

Scalar sphNFW::calcSurfaceProfile(Scalar r) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC;
	return calcScaledSurfaceProfile(r, scale);
}

void sphNFW::calcSurfaceProfile(ScalarArray1DRef r_array, ScalarArray1DRef out_array) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC;
	calcScaledSurfaceProfile(r_array, scale, out_array);
}

void sphNFW::calcSurfaceProfile(ScalarArray2DRef r_array, ScalarArray2DRef out_array) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC;
	calcScaledSurfaceProfile(r_array, scale, out_array);
}

Scalar sphNFW::calcConvergence(Scalar r, Scalar sourceSigmaC) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC/sourceSigmaC;
	return calcScaledSurfaceProfile(r, scale);
}

void sphNFW::calcConvergence(ScalarArray1DRef r_array, Scalar sourceSigmaC, ScalarArray1DRef out_array) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC/sourceSigmaC;
	calcScaledSurfaceProfile(r_array, scale, out_array);
}

void sphNFW::calcConvergence(ScalarArray2DRef r_array, Scalar sourceSigmaC, ScalarArray2DRef out_array) const
{
	Scalar scale = 2.0*rs_rhoC_deltaC/sourceSigmaC;
	calcScaledSurfaceProfile(r_array, scale, out_array);
}

Scalar sphNFW::calcShear(Scalar r, Scalar sourceSigmaC) const
{
	Scalar scale = rs_rhoC_deltaC/sourceSigmaC;
	return calcScaledShear(r, scale);
}

void sphNFW::calcShear(ScalarArray1DRef r_array, Scalar sourceSigmaC, ScalarArray1DRef out_array) const
{
	Scalar scale = rs_rhoC_deltaC/sourceSigmaC;
	calcScaledShear(r_array, scale, out_array);
}

void sphNFW::calcShear(ScalarArray2DRef r_array, Scalar sourceSigmaC, ScalarArray2DRef out_array) const
{
	Scalar scale = rs_rhoC_deltaC/sourceSigmaC;
	calcScaledShear(r_array, scale, out_array);
}
