#include <mpi.h>
#include <cmath>
#include "common.h"
#include "constants.h"
#include "mcmc.h"
#include "triaxNFW.h"

struct CatalogEntry {
	/* +00 */ double x;
	/* +08 */ double y;
	/* +16 */ double eps1;
	/* +24 */ double eps2;
	/* +32 */ double zs;
	/* +40 */ double Ds;
	/* +48 */
};

static CatalogEntry* gal_catalog;
static Vector2Array1D gal_xy(nullptr);
static ScalarArray1D gal_sigmaC(nullptr);
static int num_galaxies;
static triaxNFW nfwModel;
static Scalar zl, Dl, rhoC;
static ScalarArray1D calculated_kappas(nullptr);
static ScalarArray1D calculated_gamma1s(nullptr);
static ScalarArray1D calculated_gamma2s(nullptr);

static inline Scalar lnlikelihood(Scalar eps1, Scalar eps2, Scalar kappa, Scalar gamma1, Scalar gamma2)
{
	// TODO: calculate and return log likelihood
}

static void read_galaxy_catalog()
{
	// TODO: Allocate gal_catalog, open catalog file, and read its contents into gal_catalog and set num_galaxies
}

static void populate_global_arrays_from_galaxy_catalog()
{
	calculated_kappas.reset(new ScalarArray1Dobj(num_galaxies));
	calculated_gamma1s.reset(new ScalarArray1Dobj(num_galaxies));
	calculated_gamma2s.reset(new ScalarArray1Dobj(num_galaxies));
	gal_xy.reset(new ScalarArray1Dobj(num_galaxies));
	gal_sigmaC.reset(new Vector2Array1Dobj(num_galaxies));
	Vector2* xy = gal_xy->v;
	Scalar* p_sigmaC = gal_sigmaC->v;
	for (int n = 0; n < num_galaxies; n++, xy++, p_sigmaC++) {
		xy->x = gal_catalog[n].x;
		xy->y = gal_catalog[n].y;
		*p_sigmaC = nfwModel.calcSigmaC(gal_catalog[n].zs, gal_catalog[n].Ds);
	}
}

// The following two functions will be called from Fortran / CosmoMC

void triaxUTD_setup(double _zl, double _Dl, double _rhoC)
{
	zl=_zl;
	Dl=_Dl;
	rhoC=_rhoC;
	nfwModel.setParameters(1, 1, 1, 1, 1, 0, 0, zl, Dl, rhoC); // Set zl, Dl, and rhoC so calcSigmaC() can be used
	read_galaxy_catalog();
	populate_global_arrays_from_galaxy_catalog();
}

double triaxUTD_lnlikelihood(double c, double r200, double a, double b, double phi, double theta)
{
	nfwModel.setParameters(c, r200, NAN, a, b, theta, phi, zl, Dl, rhoC);
	nfwModel.calcConvergenceShear(gal_xy, gal_sigmaC, calculated_kappas, calculated_gamma1s, calculated_gamma2s);

	Scalar lnlk = 0.0;
	CatalogEntry* galaxy = gal_catalog;
	Scalar* p_kappa = calculated_kappas->v;
	Scalar* p_gamma1 = calculated_gamma1s->v;
	Scalar* p_gamma2 = calculated_gamma2s->v;
	for (int n = 0; n < num_galaxies; n++, galaxy++, p_kappa++, p_gamma1++, p_gamma2++) {
		lnlk += lnlikelihood(galaxy->eps1, galaxy->eps2, *p_kappa, *p_gamma1, *p_gamma2);
	}

	return lnlk;
}
