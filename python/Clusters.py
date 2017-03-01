import Cosmology,  math
import numpy as np
from scipy.integrate import quad
from scipy.misc import derivative
import matplotlib.pyplot as plt
#from pymc import *

'''
class for clusters
'''

#some constants
G = 4.5177E-48 #units of Mpc^3/solar mass/s^2
clight = 9.716E-15 #units of Mpc/s

#use this class for analytical spherical NFW functions
class sphNFW(object):

    def __init__(self, c = None, r200 = None, M200 = None, z=None, cosmoName=None):

        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmo = Cosmology.setCosmology(cosmoName)

        if z is None:
            print("No value set for the redshift. Setting z=0.25")
            self.z = 0.25
        else:
            self.z = z

        if c is None:
            print("No value set for concentration. Setting c = 4")
            self.c = 4.
        else:
            self.c = c

        if r200 is None and M200 is not None:
            self.M200 = M200
            self.r200 = self.M200Tor200()
        elif r200 is not None and M200 is None:
            self.r200 = r200
            self.M200 = self.r200ToM200()
        elif r200 is None and M200 is None:
            print ("No value set for r200 or M200. Setting r200=2 Mpc")
            self.r200 = 2.
            self.M200 = self.r200ToM200()
        else:
            raise Exception('You cannot set both M200 and r200')

    #convert r200 in Mpc to M200 in solar masses
    def r200ToM200(self):
        return 200*self.rhoC()*(4./3)*np.pi*(self.r200)**3

    #convert M200 in solar masses to r200 in Mpc
    def M200Tor200(self):
        return ((self.M200*3)/(200*self.rhoC()*4*np.pi))**(1./3)
    
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

    #return cluster scale radius in Mpc
    def rs(self):
        return self.r200 / self.c

    #Surface density profile. Input can be either 1D or 2D numpy ndarray or a float.
    def surfaceProfile(self, r=None):
        if r is None:
            raise Exception("Need to provide radius value or 1 or 2D numpy ndarray")

        rs = self.rs()
        k = 2*self.rhoC()*self.deltaC()

        if isinstance(r,np.ndarray):
            f = np.piecewise(r, [r>rs, r<rs,r==rs], [lambda r: (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5)), lambda r:(rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5))), lambda r: rs*k/3.])
        else:
            if r>rs:
                f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5))
            elif r<rs:
                f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5)))
            else:
                f = rs*k/3.
        return f


    #calculate convergence on grid or as a profile. Function of radius in Mpc
    def convergence(self, r=None, zs=None):
        if r is None:
            raise Exception("Error. Need to provide radius value or 1 or 2D numpy ndarray")
        if zs is None:
            raise Exception("Error. Need to provide source redshift (zs)")

        return self.surfaceProfile(r) / self.sigmaC(zs)

    
    #calculate shear as a profile (1D numpy array) or single value. Function of radius in Mpc.
    def shear(self, r=None, zs=None):
        rs = self.rs()
        k = self.rhoC()*self.deltaC()
        sigmaC = self.sigmaC(zs)
        x = r/rs

        if isinstance(r,np.ndarray):
            f = np.piecewise(x, [x>1., x<1., x==1.], [lambda x: (rs*k)*((8*np.arctan(((x-1)/(1+x))**0.5)/(x**2*(x**2-1)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctan(((x-1)/(1+x))**0.5)/((x**2-1)**1.5))), \
                    lambda x: (rs*k)*((8*np.arctanh(((1-x)/(1+x))**0.5)/(x**2*(1-x**2)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctanh(((1-x)/(1+x))**0.5)/((x**2-1)*(1-x**2)**0.5))), \
                    lambda x: (rs*k)*(10./3+4.*np.log(0.5))])
        else:
            if (x<1.):
                f = (rs*k)*((8*np.arctanh(((1-x)/(1+x))**0.5)/(x**2*(1-x**2)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctanh(((1-x)/(1+x))**0.5)/((x**2-1)*(1-x**2)**0.5)))
            elif (x>1.):
                f = (rs*k)*((8*np.arctan(((x-1)/(1+x))**0.5)/(x**2*(x**2-1)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctan(((x-1)/(1+x))**0.5)/((x**2-1)**1.5)))
            else:
                f = (rs*k)*(10./3+4.*np.log(0.5))

        #print f/sigmaC
        #print isinstance(f,np.ndarray)
        return f/sigmaC

    #calculate reduced shear as a profile. Function of radius in Mpc. If you're putting in x and y, they should be 1D numpy arrays of equal length
    def reducedShear(self,r=None,x=None,y=None,zs=None,tanFlag=True):
        if not tanFlag and (x is None or y is None):
            raise Exception('if you want g1 and g2 you need to give x and y coords')
        if r is not None and (x is not None or y is not None):
            raise Exception('stick to one coord system at a time')

        #print gmag
        if tanFlag:
            gtan = self.shear(r,zs) / (1-self.convergence(r,zs))
            return gtan
        if not tanFlag:
            r = (x**2.+y**2.)**0.5
            gtan = self.shear(r,zs) / (1-self.convergence(r,zs))
            locPhi = np.arctan2(y,x)
            g1 = -gtan*np.cos(2.*locPhi)
            g2 = -gtan*np.sin(2.*locPhi)
            #return np.concatenate((g1,g2),axis=0)
            return (g1,g2)


    def magnification(self,r=None,zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        kappa = self.convergence(r=r,zs=zs)
        gammaMag = self.shear(r=r,zs=zs)

        return 1/((1-kappa)**2-gammaMag**2) 

    #generate a catalog of galaxies at redshift zs with random locations and ellipticities lensed by spherical NFW cluster.
    #does not include any rejection techniques as of now. ndens is galaxy number density in arcmin^-2 and sige is dispersion of intrinsic ellip.
    #fov is the physical size of the grid over which you wanna make the catalog.

    #NEED TO ADD RETURNING G1 AND G2 FOR FITTING PURPOSES
    def galCatGen(self, zs=None, ndens=None, sige=None, fov=None,beta=0.5, tanFlag=True, randseed=None,cut = None):
        if ndens is None:
            raise Exception("No galaxy number density given.")
        if sige is None:
            raise Exception("No intrinsic ellipticity dispersion given.")
        if fov is None:
            raise Exception("No fov given.")

        #calculate the total number of galaxies in the fov
        dA = self.cosmo.angularDiameterDistance(self.z) / (self.cosmo.h)
        arcMinPerMpc = (1/dA)*10800/np.pi
        areaSqArcmin = (fov**2.)*(arcMinPerMpc**2.) #area of fov in arcmin^2
        Ngal = ndens * areaSqArcmin #gives total number of galaxies in fov

        np.random.seed(randseed)

        Ngal = int(np.random.poisson(Ngal))

        #generate random galaxy locations
        xphys = np.random.uniform(low = -fov/2.,high = fov/2., size = Ngal)
        yphys = np.random.uniform(low = -fov/2.,high = fov/2., size = Ngal)
        rphys = (xphys**2+yphys**2)**0.5

        gtan = self.reducedShear(r=rphys,zs=zs)
        #epstan = gtan+np.random.normal(0.,sige,Ngal)
        locPhi = np.arctan2(yphys,xphys)
        g1 = -gtan*np.cos(2.*locPhi)
        g2 = -gtan*np.sin(2.*locPhi)
    
        noise1 = np.random.normal(0.,sige,Ngal)/2**0.5
        noise2 = np.random.normal(0.,sige,Ngal)/2**0.5
        mu = self.magnification(r=rphys,zs=zs)

        #galaxy rejection to account for magnification of galaxy number density
        eta = np.random.rand(Ngal)
        keepVals = np.where((mu)**(beta-1)>=eta)
        #gtan = gtan[keepVals]
        rphys = rphys[keepVals]
        xphys = xphys[keepVals]
        yphys = yphys[keepVals]
        locPhi = locPhi[keepVals]
        noise1 = noise1[keepVals]
        noise2 = noise2[keepVals]
        mu = mu[keepVals]
        g1 = g1[keepVals]
        g2 = g2[keepVals]

        eps1 = (noise1+g1)/(1.+g1*noise1)
        eps2 = (noise2+g2)/(1.-g2*noise2)

        epsmag = (eps1**2.+eps2**2.)**0.5
        if cut is 'eps':
            keepVals = np.where(epsmag < 0.5)
            xphys = xphys[keepVals]
            yphys = yphys[keepVals]
            rphys = rphys[keepVals]
            eps1 = eps1[keepVals]
            eps2 = eps2[keepVals]
            locPhi = locPhi[keepVals]
        elif cut is 'mag':
            keepVals = np.where(mu < 1.8)
            xphys = xphys[keepVals]
            yphys = yphys[keepVals]
            rphys = rphys[keepVals]
            eps1 = eps1[keepVals]
            eps2 = eps2[keepVals]
            locPhi = locPhi[keepVals]
        
        Ngal = len(eps1) #readjust Ngal count

        #add noise


        if tanFlag:
            epstan = -(eps1*np.cos(2.*locPhi) + eps2*np.sin(2.*locPhi))
            return (rphys, epstan)
        else:
            return (xphys, yphys, eps1, eps2)


'''
Class for handling the simulations and doing calculations. To initialize for a cluster, you first need to read in the surface mass densities from
the .fits files and turn them into numpy arrays.

Input parameters when initializing:
sigma: 2D numpy array of surface mass density map (solar masses/Mpc^2)
z: cluster redshift (currently have only been using z=0.25)
fov: field of view of the surface mass density map. For now we've been using boxes 30 Mpc on a side, so set fov=30.
cosmoName: The cosmology used. Our sims use WMAP7-ML.
'''
class simLens(object):
    #initialize cluster. requires a surface density map, (solar masses/Mpc^2) redshift,
    # the field of view size (Mpc), and the cosmology in which it was generated.
    def __init__(self, sigma=None, z=None, fov=None, cosmoName=None):
        if sigma is None:
            raise Exception("Error. Need to include 2D surface density map (solar masses/Mpc^2) as ndarray.")
        else:
            self.sigma = sigma

        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmo = Cosmology.setCosmology(cosmoName)

        if z is None:
            print("No value set for the redshift. Setting z=0.25")
            self.z = 0.25
        else:
            self.z = z

        if fov is None:
            print("No field of view size set. Setting fov=30 Mpc")
            self.fov = 30.
        else:
            self.fov = fov

    #return critical surface density for cluster and source redshifts in solar masses/Mpc^2
    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        zl = self.z
        Dl = self.cosmo.angularDiameterDistance(zl) / self.cosmo.h
        Ds = self.cosmo.angularDiameterDistance(zs) / self.cosmo.h
        Dls = (1/(1+zs))*(Ds*(1+zs)-Dl*(1+zl))
        return clight**2*Ds/(4*math.pi*G*Dls*Dl)

    #returns convergence map on grid
    def convergence(self,zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
            
        return self.sigma/self.sigmaC(zs)

    #returns shear map grid calculated with FFT. Set polar=True to get (shear magnitude, shear orientation).
    #Otherwise returns (gamma1,gamma2).
    def shear(self,zs=None, polar=False):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        #fft on convergence grid
        kappaF = np.fft.fft2(self.convergence(zs))
        pix = self.fov/kappaF.shape[0]

        #get the frequencies using pixel size as spacing
        freqx = np.fft.fftfreq(kappaF.shape[0],d=pix)
        freqy = np.fft.fftfreq(kappaF.shape[1],d=pix)
        freqY, freqX = np.meshgrid(freqx, freqy)

        #initialize and then calculate shear in fourier space
        gamma1F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)
        gamma2F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)

        gamma1F = kappaF*(freqX**2-freqY**2)/(freqX**2+freqY**2)
        gamma2F = 2.*kappaF*(freqX*freqY)/(freqX**2+freqY**2)

        #replace bad elements
        gamma1F = np.nan_to_num(gamma1F)
        gamma2F = np.nan_to_num(gamma2F)

        #fft back to position space
        gamma1 = np.fft.ifft2(gamma1F).real
        gamma2 = np.fft.ifft2(gamma2F).real

        gammaMag = (gamma1**2+gamma2**2)**(1./2)
        gammaPhi = (1./2)*np.arctan2(gamma2,gamma1)
        #gammaPhi = np.pi/2. - gammaPhi
        #gamma1 = gammaMag*np.cos(2.*gammaPhi)
        #gamma2 = gammaMag*np.sin(2.*gammaPhi)
        #print gammaMag
        #print gamma2
        #print gammaPhi

        if polar:
            return (gammaMag, gammaPhi)
        else:
            return (gamma1, gamma2)


    def magnification(self,zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        kappa = self.convergence(zs=zs)
        gammaMag, gammaPhi = self.shear(zs=zs,polar=True)

        return 1/((1-kappa)**2-gammaMag**2) 

    #calculate reduced shear on a grid for sim clusters. Set polar=True to get (g magnitude, g orientation).
    #Otherwise returns (g1,g2)
    def reducedShear(self,zs=None, polar=False, tan=False, retkappa=False):
        if zs is None:
            raise Exception("Error. Need to provide source redshift (zs)")
        if polar and tan:
            raise Eception("Error. Can't set both polar and tan as True.")

        gamma1,gamma2 = self.shear(zs)
        kappa = self.convergence(zs)
        g1 = gamma1/(1-kappa)
        g2 = gamma2/(1-kappa)
        npix = g1.shape[0]

        gMag = (g1**2 + g2**2)**0.5
        gPhi = 0.5*np.arctan2(g2,g1)

        if polar and not retkappa:
            return (gMag,gPhi)
        elif tan:
            xphys = np.linspace(-self.fov/2.,self.fov/2.,npix)
            yphys = np.linspace(-self.fov/2.,self.fov/2.,npix)
            yyphys,xxphys = np.meshgrid(xphys,yphys)
            locPhi = np.arctan2(yyphys,xxphys)
            gtan = -(g1*np.cos(2.*locPhi) + g2*np.sin(2.*locPhi))
            return (xxphys,yyphys,gtan)
        elif not tan and not retkappa:
            return (g1,g2)
        elif polar and retkappa:
            return (kappa,gMag,gPhi)

    '''
    generate a catalog of galaxies at redshift zs with random locations and ellipticities lensed by cosmo-OWLS cluster.
    does not include any rejection techniques as of now.

    zs: source redshift. We normally use zs=1.
    ndens: galaxy number density in arcmin^-2. Usually set to 30.
    sige: dispersion of intrinsic ellip. Usually set to 0.25 or so, but you can set it to something really small for debugging and testing

    #set tanFlag to true to return only tangential component of the shear
    '''
    def galCatGen(self, zs=None, ndens=None, sige=None, randseed = None, beta=0.5, tanFlag=False, gMag=None, gPhi=None, mu=None,retkappa=False, cut=None):
        if ndens is None:
            raise Exception("Error. No galaxy number density given.")
        if sige is None:
            raise Exception("Error. No intrinsic ellipticity dispersion given.")

        #calculate the total number of galaxies in the fov
        dA = self.cosmo.angularDiameterDistance(self.z) / (self.cosmo.h)
        arcMinPerMpc = (1/dA)*10800/np.pi
        areaSqArcmin = (self.fov**2.)*(arcMinPerMpc**2.) #area of fov in arcmin^2
        Ngal = ndens * areaSqArcmin #gives total number of galaxies in fov

        #calculate reduced shear everywhere. Have to adjust gPhi for fitting bc python is weird about how it stores arrays.
        if gMag is None and not retkappa:
            gMag,gPhi = self.reducedShear(zs=zs, polar=True)
            #g1,g2 = self.reducedShear(zs=zs)
            mu = self.magnification(zs=zs)
            #gPhi = np.pi/2. - gPhi
        elif gMag is None and retkappa:
            kappa,gMag,gPhi = self.reducedShear(zs=zs,polar=True,retkappa=True)
            mu = self.magnification(zs=zs)

        g1 = gMag*np.cos(2.*gPhi)
        g2 = gMag*np.sin(2.*gPhi)
        pix = self.fov/gMag.shape[0]
        npix = gMag.shape[0]

        #set seed if requested for constant galaxy catalogs
        np.random.seed(randseed)

        #poisson adjustment to Ngal
        Ngal = int(np.random.poisson(Ngal))

        #generate random galaxy locations
        xind = np.random.randint(low = 0,high = npix, size = Ngal)
        yind = np.random.randint(low = 0, high = npix, size = Ngal)

        #convert pixels to physical locations
        xphys = (xind*pix) - (self.fov/2.)
        yphys = (yind*pix) - (self.fov/2.)
        rphys = (xphys**2 + yphys**2)**0.5

        #get the subset of g1, g2, and mag for the galaxy locations
        g1sub = g1[xind,yind]
        g2sub = g2[xind,yind]
        musub = mu[xind,yind]
        noise1 = np.random.normal(0.,sige,Ngal)/(2**0.5)
        noise2 = np.random.normal(0.,sige,Ngal)/(2**0.5)
        if retkappa:
            kappasub = kappa[xind,yind]

        eps1 = (noise1+g1sub)/(1.+g1sub*noise1)
        eps2 = (noise2+g2sub)/(1.-g2sub*noise2)
        epsmag = (eps1**2.+eps2**2.)**0.5

        #write the galcat out before mag rejection
        arrayToWrite = [xphys, yphys, eps1,eps2]
        np.savetxt('galcat_lownoise_1000_x.dat',np.transpose(arrayToWrite),delimiter=' ',fmt = '%.5f')

        #galaxy rejection to account for magnification of galaxy number density
        eta = np.random.rand(Ngal)
        keepVals = np.where((musub)**(beta-1)>=eta)
        eps1 = eps1[keepVals]
        eps2 = eps2[keepVals]
        epsmag = epsmag[keepVals]
        xphys = xphys[keepVals]
        yphys = yphys[keepVals]
        rphys = rphys[keepVals]
        musub = musub[keepVals]
        if retkappa:
            kappasub = kappasub[keepVals]


        #add noise
        #eps1 = g1sub+noise1
        #eps2 = g2sub+noise2
        
        if cut is 'eps':
            print 'eps cut'
            keepVals = np.where(epsmag < 0.5)
            epsmag = epsmag[keepVals]
            xphys = xphys[keepVals]
            yphys = yphys[keepVals]
            rphys = rphys[keepVals]
            eps1 = eps1[keepVals]
            eps2 = eps2[keepVals]
        elif cut is 'mag':
            print 'mag cut'
            keepVals = np.where(musub < 1.8)
            epsmag = epsmag[keepVals]
            xphys = xphys[keepVals]
            yphys = yphys[keepVals]
            rphys = rphys[keepVals]
            eps1 = eps1[keepVals]
            eps2 = eps2[keepVals]

        Ngal = len(g1sub) #readjust Ngal count

        #arrayToWrite = [xphys,yphys,kappasub,eps1,eps2]
        #np.savetxt('galcat.dat',np.transpose(arrayToWrite),delimiter=' ',fmt='%.5f')

        #eps1 = g1sub
        #eps2 = g2sub

        #fig1 = plt.figure()
        #plt.hist(noise1,bins=30)
        #plt.hist(eps1,bins=30)
        #fig2 = plt.figure()
        #plt.hist(noise2,bins=30)
        #fig3 = plt.figure()
        #plt.show()

        #print g1sub.shape
        #print xphys
        #print yphys

        if tanFlag:
            locPhi = np.arctan2(yphys,xphys)
            #print locPhi
            epstan = -(eps1*np.cos(2.*locPhi) + eps2*np.sin(2.*locPhi))
            #epstan = (eps1**2.+eps2**2.)**0.5
            #epstan = epstan + noise1
            return (rphys, epstan)
        else:
            return (xphys, yphys, eps1, eps2)
'''
the purpose of this class is to perform the numerical and analytical calculations for the triaxial NFW model for the purposes of generating ideal triaxial NFW clusters and to provide the model for fitting.

currently this does not have the galaxy catalog generator as in the other two classes but this is to be added.

to initialize, you need the 6 main parameters (c, r200 or M200, the axis ratios a and b, and the orientation angles theta and phi) as well as the lens redshift z and the name of the cosmology you're using

see feroz and hobson 2012 for some mathematical details of model
'''
class triaxNFW(object):

    def __init__(self, c = None, r200 = None, M200 = None, a=None, b = None, theta=None, phi=None, z=None, cosmoName=None):

        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmoName = cosmoName
        self.cosmo = Cosmology.setCosmology(cosmoName)

        if z is None:
            print("No value set for the redshift. Setting z=0.25")
            self.z = 0.25
        else:
            self.z = z

        if c is None:
            print("No value set for concentration. Setting c = 4")
            self.c = 4.
        else:
            self.c = c

        if theta is None:
            print("No value set for theta. Setting theta=0")
            self.theta = 0.
        else:
            self.theta = theta

        if phi is None:
            print("No value set for phi. Setting phi=0")
            self.phi = 0.
        else:
            self.phi = phi

        if a is None or b is None:
            print("No value set for a or b or for both. Setting a=b=1 (spherical case)")
            self.a = 1.
            self.b = 1.
        else:
            self.a = a
            self.b = b

        if r200 is None and M200 is not None:
            self.M200 = M200
            self.r200 = self.M200Tor200()
        elif r200 is not None and M200 is None:
            self.r200 = r200
            self.M200 = self.r200ToM200()
        elif r200 is None and M200 is None:
            print ("No value set for r200 or M200. Setting r200=2 Mpc")
            self.r200 = 2.
            self.M200 = self.r200ToM200()
        else:
            raise Exception('You cannot set both M200 and r200')

        #make a spherical NFW cluster for later use
        self.nfwcluster = sphNFW(c=self.c,r200=self.r200,z=self.z,cosmoName=self.cosmoName)

    #convert r200 in Mpc to M200 in solar masses. MAKE SURE CONVERSIONS ARE RIGHT BC OF SHAPE
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

    #critical surface density for getting convergence
    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        zl = self.z
        Dl = self.cosmo.angularDiameterDistance(zl) / self.cosmo.h
        Ds = self.cosmo.angularDiameterDistance(zs) / self.cosmo.h
        Dls = (1/(1+zs))*(Ds*(1+zs)-Dl*(1+zl))
        return clight**2*Ds/(4*math.pi*G*Dls*Dl)

    #some of the constants that need to be calculated
    def returnfABC(self):
        f = (np.sin(self.theta)**2)*((np.cos(self.phi)**2)/(self.a**2) + (np.sin(self.phi)**2)/(self.b**2)) + np.cos(self.theta)**2
        A = (np.cos(self.theta)**2)*((np.sin(self.phi)**2)/(self.a**2) + (np.cos(self.phi)**2)/(self.b**2)) + (np.sin(self.theta)**2)/((self.a**2)*(self.b**2))
        B = np.cos(self.theta)*np.sin(2*self.phi)*((1/self.a)**2 - (1/self.b)**2)
        C = (np.sin(self.phi)**2)/(self.b**2) + (np.cos(self.phi)**2)/(self.a**2)

        return (f, A, B, C)

    #calculate projected axis ratios
    def returnq(self, f=None, A=None, B=None, C=None):
        if f is None or A is None or B is None or C is None:
            f,A,B,C = self.returnfABC()

        qX = (2*f/(A+C - ((A-C)**2+B**2)**.5))**.5
        qY = (2*f/(A+C + ((A-C)**2+B**2)**.5))**.5

        #axis ratio of elliptical contours
        q = qY/qX

        return (qX, qY, q)

    #another angle for use later when rotating back to original coordinate system
    def returnpsi(self, A=None,B=None,C=None):
        if A is None or B is None or C is None:
            f, A, B, C = self.returnfABC()

        return np.arctan2(B, A-C)

    #triaxial radius squared
    def zetasq(self, u, x, y, q=None, prescaled=False):
        if prescaled and q is None:
            qx,qy,q = self.returnq()
        elif prescaled is False:
            qx,qy,q = self.returnq()
            x=x/qx
            y=y/qx

        return u*(x**2 + (y**2)/(1-(1-q**2)*u))

    #functions for calculating the necessary integrals
    def Jintegrand(self,u,x,y,q,f,nfwcluster, zs, n):
        zeta = (self.zetasq(u,x,y,q,prescaled=True))**0.5
        return (1/f**0.5)*nfwcluster.convergence(zeta,zs)/(1-(1-q**2)*u)**(n+0.5)


    def Kintegrand(self, u,x,y,q,f,nfwcluster,zs,n):
        zeta = (self.zetasq(u,x,y,q,prescaled=True))**0.5
        delkappa = derivative(nfwcluster.convergence,zeta,dx=1e-5,args=(zs,))
        return (1/f**0.5)*u*delkappa/(zeta*(1-(1-q**2)*u)**(n+0.5))

    def JKintegrals(self, x, y, q, f, zs): #these x and y are assumed to be prescaled by qx
        J0, errJ0 = quad(self.Jintegrand,0,1, args=(x,y,q,f,self.nfwcluster,zs,0))
        J1, errJ1 = quad(self.Jintegrand,0,1, args=(x,y,q,f,self.nfwcluster,zs,1))
        K0, errK0 = quad(self.Kintegrand,0,1, args=(x,y,q,f,self.nfwcluster,zs,0))
        K1, errK1 = quad(self.Kintegrand,0,1, args=(x,y,q,f,self.nfwcluster,zs,1))
        K2, errK2 = quad(self.Kintegrand,0,1, args=(x,y,q,f,self.nfwcluster,zs,2))
        return (J0,J1,K0,K1,K2)

    #calculate the second derivatives of the potential
    def potential(self,x,y,q,f,zs): #these x and y are assumed to be prescaled by qx
        J0,J1,K0,K1,K2 = self.JKintegrals(x,y,q,f,zs)
        potxx = q*x**2*K0 + q*J0
        potyy = q*y**2*K2 + q*J1
        potxy = q*x*y*K1
        return (potxx,potyy,potxy)

    #calculate the convergence and shear in intermediate coordinate system
    def transconvshear(self,x,y,q,f,zs): #convergence and shear calculated before switch to original coords. x and y are assumed to be prescaled by qx
        potxx,potyy,potxy = self.potential(x,y,q,f,zs)
        kappat = 0.5*(potxx+potyy)
        gamma1t = 0.5*(potxx-potyy)
        gamma2t = potxy

        return (kappat,gamma1t,gamma2t)

    #return back to intermediate coordinate system for convergence and shear. this may need to be adjusted to loop over list of (x,y) instead of looping over all y for every x coords
    def convshear(self,x,y,zs): #the input x and y are NOT SCALED BY QX, THESE ARE THE ORIGINAL COORDINATES YOU WANT THE MAPS IN
        f,A,B,C = self.returnfABC()
        q,qx,qy = self.returnq(f,A,B,C)

        #scaling x and y
        x = x/qx
        y = y/qx

        x = x.tolist()
        y = y.tolist()
        xnpix = len(x)
        ynpix = len(y)

        kappa = []
        gamma1t = []
        gamma2t = []
        for ix,xval in enumerate(x):
            for iy,yval in enumerate(y):
                kappatemp, gamma1temp, gamma2temp = self.transconvshear(xval,yval,q,f,zs)
                kappa.append(kappatemp)
                gamma1t.append(gamma1temp)
                gamma2t.append(gamma2temp)

        kappa = np.asarray(kappa).reshape([xnpix,ynpix])
        gamma1t = np.asarray(gamma1t).reshape([xnpix,ynpix])
        gamma2t = np.asarray(gamma2t).reshape([xnpix,ynpix])
        alpha = np.arctan2(gamma2t,gamma1t)
        psi = self.returnpsi(A,B,C)
        gammaMag = (gamma1t**2+gamma2t**2)**0.5
        gamma1 = gammaMag*np.cos(2*(alpha+psi))
        gamma2 = gammaMag*np.sin(2*(alpha+psi))

        return kappa, gamma1, gamma2
