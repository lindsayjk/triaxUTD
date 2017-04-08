#pragma once
#include <cmath>
#include "common.h"
#include "constants.h";
#include "Cosmology.h"

// Base class for deriving various lensing models.
class LensModel {
private:
	LensModel() {}

public:
	// Either r200 or M200 may be NAN but not both.
	LensModel(Scalar c, Scalar r200, Scalar M200, Scalar z, const Cosmology* cosmo);
	virtual ~LensModel();

	// convert r200 in Mpc to M200 in solar masses
	Scalar r200ToM200() const
	{
		return 200*get_rhoC()*(4./3)*M_PI*pow(r200, 3.0);
	}

	// convert M200 in solar masses to r200 in Mpc
	Scalar M200tor200() const
	{
		return cbrt((M200*3)/(200*get_rhoC()*4*M_PI));
	}

	// return cluster scale radius in Mpc
	Scalar get_rs() const
	{
		return r200 / c;
	}

	// return critical density at cluster redshift in units of solar masses/Mpc^3
	Scalar get_rhoC() const
	{
		return cosmo->rho_c(z)*(1.E9)*cosmo->h2();
	}

	// return delta_C for cluster concentration
	Scalar get_deltaC() const
	{
		return (200./3.)*pow(c,3.0)/(log(1.0+c)-(c/(1.0+c)));
	}

	// return critical surface density for cluster and source redshifts in solar masses/Mpc^2
	Scalar get_sigmaC(Scalar z_source) const
	{
		Scalar z_lens = z;
		Scalar Dl = cosmo->angularDiameterDistance(z_lens) / cosmo->h();
		Scalar Ds = cosmo->angularDiameterDistance(z_source) / cosmo->h();
		Scalar Dls = (1.0/(1.0+zs))*(Ds*(1.0+zs)-Dl*(1.0+zl));
		return Constants.clight2*Ds/(4*M_PI*Constants.G*Dls*Dl);
	}

protected:
	const Cosmology* cosmo;
	Scalar c, r200, M200, z;
};
