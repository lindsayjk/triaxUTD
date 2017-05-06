#pragma once

struct CatalogEntry {
	/* +00 */ double x;
	/* +08 */ double y;
	/* +16 */ double eps1;
	/* +24 */ double eps2;
	/* +32 */ double zs;
	/* +40 */ double Ds;
	/* +48 */
};

extern "C" {
	void triaxUTD_setup();
	double triaxUTD_lnlikelihood(double c, double r200, double a, double b, double phi, double theta);
}
