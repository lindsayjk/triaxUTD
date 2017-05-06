#pragma once

extern "C" {
	void triaxUTD_setup();
	double triaxUTD_lnlikelihood(double c, double r200, double a, double b, double phi, double theta);
}
