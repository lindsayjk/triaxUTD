#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstring>
#include "common.h"
#include "constants.h"
#include "mcmc.h"
#include "mpi_logging.h"
#include "triaxNFW.h"

using namespace std;

static string gal_catalog_path;
static CatalogEntry* gal_catalog = nullptr;
static Vector2Array1D gal_xy(nullptr);
static Vector2Array1D gal_eps(nullptr);
static ScalarArray1D gal_sigmaC(nullptr);
static int num_galaxies;
static triaxNFW nfwModel;
static Scalar zl, Dl, rhoC;
static ScalarArray1D calculated_kappas(nullptr);
static ScalarArray1D calculated_gamma1s(nullptr);
static ScalarArray1D calculated_gamma2s(nullptr);

static inline Scalar lnlikelihood(Scalar eps1, Scalar eps2, Scalar kappa, Scalar gamma1, Scalar gamma2)
{
	const double sige= 0.25;
	const double sige2 = 0.0625;

	double g1model = gamma1/(1-kappa);
	double g2model = gamma2/(1-kappa);

	complex<double> gmod(g1model,g2model);
	double gabs2 = pow(abs(gmod),2);
	complex<double> gmodconj = conj(gmod);
	complex<double> epsdat(eps1, eps2);


	complex<double> epssrc = (epsdat - gmod)/(1-gmodconj*epsdat);
	double epssrcmag = abs(epssrc);
	double fac1 = -pow(epssrcmag,2)/sige2; //doesn't need to be ln'd, other 3 terms do
	double fac2 = pow((gabs2-1),2);
	double fac3 = M_PI*sige2*(1-exp(-1/sige2));
	double fac4 = pow(abs(epsdat*gmodconj - 1),4);

	double loglike = fac1 + log(fac2) - log(fac3) - log(fac4);

	return loglike;

}

static void read_galaxy_catalog()
{
	size_t fsize;
	ifstream f(gal_catalog_path, ifstream::binary | ifstream::ate);
	if(f.fail()) throw_line("Galaxy catalog file could not be opened!");
	fsize = f.tellg();
	f.seekg(0);
	num_galaxies = fsize / sizeof(CatalogEntry);
	if(gal_catalog) delete gal_catalog;
	gal_catalog = new CatalogEntry[num_galaxies];
	f.read((char*)gal_catalog, fsize);
	f.close();
}

static void populate_global_arrays_from_galaxy_catalog()
{
	const Scalar cutoff_radius = 0.5;
	int num_valid_galaxies = 0;
	int n;
	for(n = 0;n < num_galaxies; n++) {
		Scalar x = gal_catalog[n].x;
		Scalar y = gal_catalog[n].y;
		if(sqrt((x*x)+(y*y))>=cutoff_radius) {
			num_valid_galaxies++;
		} else {
			gal_catalog[n].x=NAN;
			gal_catalog[n].y=NAN;
		}
	}
	mpi_log(nullptr, "Filtered out %d/%d galaxies within %g Mpc", num_galaxies-num_valid_galaxies, num_galaxies, cutoff_radius);
	calculated_kappas.reset(new ScalarArray1Dobj(num_valid_galaxies));
	calculated_gamma1s.reset(new ScalarArray1Dobj(num_valid_galaxies));
	calculated_gamma2s.reset(new ScalarArray1Dobj(num_valid_galaxies));
	gal_xy.reset(new Vector2Array1Dobj(num_valid_galaxies));
	gal_eps.reset(new Vector2Array1Dobj(num_valid_galaxies));
	gal_sigmaC.reset(new ScalarArray1Dobj(num_valid_galaxies));
	Vector2* xy = gal_xy->v;
	Vector2* eps = gal_eps->v;
	Scalar* p_sigmaC = gal_sigmaC->v;
	for (n = 0; n < num_galaxies; n++) {
		if(std::isnan(gal_catalog[n].x)) continue;
		xy->x = gal_catalog[n].x;
		xy->y = gal_catalog[n].y;
		eps->x = gal_catalog[n].eps1;
		eps->y = gal_catalog[n].eps2;
		*p_sigmaC = nfwModel.calcSigmaC(gal_catalog[n].zs, gal_catalog[n].Ds);
		xy++;
		eps++;
		p_sigmaC++;
	}
	num_galaxies = num_valid_galaxies;

	xy = gal_xy->v;
	eps = gal_eps->v;
	p_sigmaC = gal_sigmaC->v;
	for(n = 0;n < num_galaxies;n++) {
		if(sqrt(xy->x*xy->x+xy->y*xy->y)<cutoff_radius) {
			mpi_log(nullptr, "Found galaxy %d/%d inside cutoff radius after filtering! x=%g, y=%g, eps=(%g,%g) sigmaC=%g",n+1,num_galaxies,xy->x,xy->y,eps->x,eps->y,*p_sigmaC);
		}
		xy++;
		eps++;
		p_sigmaC++;
	}
}

static void triaxUTD_terminate()
{
	static bool rethrown = false;
	try {
		if (!rethrown) {
			rethrown = true;
			throw;
		}
	}
	catch (const std::exception& e) {
		cerr << "C++ exception: " << e.what() << endl;
	}
	catch (...) {
		cerr << "Unknown C++ exception caught by triaxUTD_terminate" << endl;
	}
	std::abort();
}

// The following two functions will be called from Fortran / CosmoMC

void triaxUTD_setup()
{
	init_mpi_logging("triaxutd.log", "triaxutd.log");
	std::set_terminate(triaxUTD_terminate);
	init_gsl_error_handling();

	mpi_log(nullptr, "Process started");

	ifstream params_file("triaxutd.params");
	string line;
	double d;

	if(params_file.fail()) throw_line("triaxutd.params could not be opened!");

	getline(params_file, line); // First line is a comment
	getline(params_file, line); // Second line is a comment
	getline(params_file, gal_catalog_path); // Path to galaxy catalog file

	getline(params_file, line); // zl
	sscanf(line.c_str(), "%lf", &d);
	zl=d;

	getline(params_file, line); // Dl
	sscanf(line.c_str(), "%lf", &d);
	Dl=d;

	getline(params_file, line); // rhoC
	sscanf(line.c_str(), "%lf", &d);
	rhoC=d;

	mpi_log(nullptr, "triaxUTD_setup gal_catalog_path=\"%s\" zl=%f Dl=%f rhoC=%f", gal_catalog_path.c_str(), zl, Dl, rhoC);

	params_file.close();

	nfwModel.setParameters(1, 1, 1, 1, 1, 0, 0, zl, Dl, rhoC); // Set zl, Dl, and rhoC so calcSigmaC() can be used
	read_galaxy_catalog();
	populate_global_arrays_from_galaxy_catalog();
}

double triaxUTD_lnlikelihood(double c, double r200, double a, double b, double phi, double theta)
{
/*
	static int count = 0;
	count++;
	if((count%100)==0) {
		mpi_log(nullptr, "triaxUTD_lnlikelihood called with c=%f r200=%f a=%f b=%f phi=%f theta=%f\n", c, r200, a, b, phi, theta);
	}*/

	if(a>b) return /*-INFINITY*/1.0E30;
	nfwModel.setParameters(c, r200, NAN, a, b, theta, phi, zl, Dl, rhoC);
	nfwModel.calcConvergenceShear(gal_xy, gal_sigmaC, calculated_kappas, calculated_gamma1s, calculated_gamma2s);

	Scalar lnlk = 0.0;
	//CatalogEntry* galaxy = gal_catalog;
	Vector2* p_eps = gal_eps->v;
	Scalar* p_kappa = calculated_kappas->v;
	Scalar* p_gamma1 = calculated_gamma1s->v;
	Scalar* p_gamma2 = calculated_gamma2s->v;
	for (int n = 0; n < num_galaxies; n++, /*galaxy++,*/ p_eps++,  p_kappa++, p_gamma1++, p_gamma2++) {
//		lnlk += lnlikelihood(galaxy->eps1, galaxy->eps2, *p_kappa, *p_gamma1, *p_gamma2);
		lnlk += lnlikelihood(p_eps->x, p_eps->y, *p_kappa, *p_gamma1, *p_gamma2);
	}

	return -lnlk;
}

