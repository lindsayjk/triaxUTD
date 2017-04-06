#pragma once
#include <cmath>
#include "common.h"
#include "constants.h";
#include "Cosmology.h"

/*
sphNFW
	def __init__(self, c = None, r200 = None, M200 = None, z=None, cosmoName=None):

	#convert r200 in Mpc to M200 in solar masses
	def r200ToM200(self):
	return 200*self.rhoC()*(4./3)*np.pi*(self.r200)**3

	#convert M200 in solar masses to r200 in Mpc
	def M200Tor200(self):
	return ((self.M200*3)/(200*self.rhoC()*4*np.pi))**(1./3)

	#return cluster scale radius in Mpc
	def rs(self):
	return self.r200 / self.c

	#return critical density at cluster redshift in units of solar masses/Mpc^3
	def rhoC(self):
	return self.cosmo.rho_c(self.z)*1000**3*self.cosmo.h**2

	#return delta_C for cluster concentration
	def deltaC(self):
	return (200./3.)*(self.c**3)/(np.log(1+self.c)-(self.c/(1+self.c)))

	#return critical surface density for cluster and source redshifts in solar masses/Mpc^2
	def sigmaC(self, zs=None):
	if zs is None:
	raise Exception("Error. Need to set source redshift (zs)")
	zl = self.z
	Dl = self.cosmo.angularDiameterDistance(zl) / self.cosmo.h
	Ds = self.cosmo.angularDiameterDistance(zs) / self.cosmo.h
	Dls = (1/(1+zs))*(Ds*(1+zs)-Dl*(1+zl))
	return clight**2*Ds/(4*math.pi*G*Dls*Dl)

	#Surface density profile. Input can be either 1D or 2D numpy ndarray or a float.
	def surfaceProfile(self, r=None):
	
	#calculate convergence on grid or as a profile. Function of radius in Mpc
	def convergence(self, r=None, zs=None):
	
	#calculate shear as a profile (1D numpy array) or single value. Function of radius in Mpc.
	def shear(self, r=None, zs=None):
	
triaxNFW
	def __init__(self, c = None, r200 = None, M200 = None, a=None, b = None, theta=None, phi=None, z=None, cosmoName=None):

	#some of the constants that need to be calculated
	def returnfABC(self):
	
	#calculate projected axis ratios
	def returnq(self, f=None, A=None, B=None, C=None):
	
	#another angle for use later when rotating back to original coordinate system
	def returnpsi(self, A=None,B=None,C=None):

	#triaxial radius squared
	def zetasq(self, u, x, y, q=None, prescaled=False):
	
	#functions for calculating the necessary integrals
	def Jintegrand(self,u,x,y,q,f,nfwcluster, zs, n):
	def Kintegrand(self, u,x,y,q,f,nfwcluster,zs,n):
	def JKintegrals(self, x, y, q, f, zs): #these x and y are assumed to be prescaled by qx

	#calculate the second derivatives of the potential
	def potential(self,x,y,q,f,zs): #these x and y are assumed to be prescaled by qx

	#calculate the convergence and shear in intermediate coordinate system
	def transconvshear(self,x,y,q,f,zs): #convergence and shear calculated before switch to original coords. x and y are assumed to be prescaled by qx

	#return back to intermediate coordinate system for convergence and shear. this may need to be adjusted to loop over list of (x,y) instead of looping over all y for every x coords
	def convshear(self,x,y,zs): #the input x and y are NOT SCALED BY QX, THESE ARE THE ORIGINAL COORDINATES YOU WANT THE MAPS IN
	*/

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
