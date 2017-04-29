#pragma once

extern "C" {
	void triaxUTD_setup(double zl, double Dl, double rhoC);
	double triaxUTD_lnlikelihood(double c, double r200, double a, double b, double phi, double theta);
}
