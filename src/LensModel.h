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
		return 200*rhoC*(4./3)*M_PI*pow(r200, 3.0);
	}

	// convert M200 in solar masses to r200 in Mpc
	Scalar M200tor200() const
	{
		return cbrt((M200*3)/(200*rhoC*4*M_PI));
	}

	// return cluster scale radius in Mpc
	Scalar get_rs() const
	{
		return rs;
	}

	// return critical density at cluster redshift in units of solar masses/Mpc^3
	Scalar get_rhoC() const
	{
		return rhoC;
	}

	// return deltaC for cluster concentration
	Scalar get_deltaC() const
	{
		return deltaC;
	}

	// return critical surface density for cluster and source redshifts in solar masses/Mpc^2
	// Performance Note: This calls into Cosmology.
	Scalar calcSigmaC(Scalar z_source) const
	{
		Scalar Ds = cosmo->angularDiameterDistance(z_source) / cosmo->h();
		Scalar Dls = (1.0/(1.0+z_source))*(Ds*(1.0+z_source)-Dl*(1.0+z));
		return Constants.clight2*Ds/(4*M_PI*Constants.G*Dls*Dl);
	}

protected:
	const Cosmology* cosmo;
	Scalar c, r200, M200, z;
	Scalar Dl;
	Scalar rs; // cluster scale radius in Mpc
	Scalar rhoC; // critical density at cluster redshift in units of solar masses/Mpc^3
	Scalar deltaC; // deltaC for cluster concentration
	Scalar rs_rhoC_deltaC; // rs, rhoC, and deltaC multiplied together
};
