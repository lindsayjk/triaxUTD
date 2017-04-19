#pragma once
#include "LensModel.h"

class sphNFW : public LensModel {
private:
	sphNFW() {}
public:
	// Either r200 or M200 may be NAN but not both.
	sphNFW(Scalar c, Scalar r200, Scalar M200, Scalar z, const Cosmology* cosmo);
	~sphNFW();

	Scalar calcSurfaceProfile(Scalar r) const;
	void calcSurfaceProfile(ScalarArray1DRef r, ScalarArray1DRef out) const;
	void calcSurfaceProfile(ScalarArray2DRef r, ScalarArray2DRef out) const;
	
	// To get sourceSigmaC, call calcSigmaC with the source redshift (calcSigmaC is defined in LensModel)

	Scalar calcConvergence(Scalar r, Scalar sourceSigmaC) const;
	void calcConvergence(ScalarArray1DRef r, Scalar sourceSigmaC, ScalarArray1DRef out) const;
	void calcConvergence(ScalarArray2DRef r, Scalar sourceSigmaC, ScalarArray2DRef out) const;

	Scalar calcShear(Scalar r, Scalar sourceSigmaC) const;
	void calcShear(ScalarArray1DRef r, Scalar sourceSigmaC, ScalarArray1DRef out) const;
	void calcShear(ScalarArray2DRef r, Scalar sourceSigmaC, ScalarArray2DRef out) const;

protected:
	Scalar calcScaledSurfaceProfile(Scalar r, Scalar scale) const;
	void calcScaledSurfaceProfile(ScalarArray1DRef r, Scalar scale, ScalarArray1DRef out) const;
	void calcScaledSurfaceProfile(ScalarArray2DRef r, Scalar scale, ScalarArray2DRef out) const;

	Scalar calcScaledShear(Scalar r, Scalar scale) const;
	void calcScaledShear(ScalarArray1DRef r, Scalar scale, ScalarArray1DRef out) const;
	void calcScaledShear(ScalarArray2DRef r, Scalar scale, ScalarArray2DRef out) const;
};
