#!/usr/bin/env python

'''
Profile curves for a ring and some spirals in the xy plane.


#Spiral Archimedean Polyskelion
#Ref : https://math.stackexchange.com/questions/651772/parametric-equations-and-specifications-of-a-triskelion-triple-spiral

'''

import math
import numpy as np
from scipy.special import fresnel, hyp2f1
from scipy.optimize import root_scalar
from scipy.integrate import quadrature


class SpiralArchimedes(object):
    '''
    Polar equation: r = b*theta

    '''
    def __init__(self, b, rmin):
        if b <= 0:
            raise ValueError('b (= %g) must be > 0.'%b)
        if rmin < 0:
            raise ValueError('rmin (= %g) must be >= 0.'%rmin)
        self.b = b
        self.rmin = rmin
        self.tmin = self.r_to_theta(self.rmin)


    def r_to_theta(self, r):
        r_ = np.asarray(r)
        if np.any(r_ < self.rmin):
            raise ValueError('r must be >= rmin (= %g).'%self.rmin)
        return r_/self.b


    def theta_to_r(self, theta):
        theta_ = np.asarray(theta)
        if np.any(theta_ < self.tmin):
            raise ValueError('theta must be >= tmin (= %g).'%self.tmin)
        return self.b*theta_


    def arclength_to_theta(self, s):
        '''
        Returns the theta parameter at arclength s. Note that here the 
        arclength s is measured from theta = tmin, where tmin is the starting
        value of theta for the spiral.

        '''
        if np.isscalar(s):
            s_ = np.array([s]); s_isscalar = True
        else:
            s_ = np.asarray(s); s_isscalar = False
        if np.any(s_ < 0):
            raise ValueError('s must be >= 0.')
        out = np.empty_like(s_)
        n = s_.size
        for i in range(n):
            sol = root_scalar(self.__class__._al2t_func, 
                    args=(self.b, self.tmin, s_[i]), 
                    method='newton', fprime=True, x0=self.tmin, xtol=1.e-8,
                    rtol=1.e-10)
            #Check convergence
            if not sol.converged:
                print('theta calculation not converged for element %d'%i)
            out[i] = sol.root
        if s_isscalar:
            return out[0]
        else:
            return out


    @staticmethod
    def _al2t_func(t, *args):
        b = args[0]; tmin = args[1]; s = args[2] 
        C = tmin*np.sqrt(1+tmin*tmin) + np.arcsinh(tmin)
        f = t*np.sqrt(1+t*t) + np.arcsinh(t) - (2*s/b + C)
        fder = 2*np.sqrt(1+t*t)
        return f, fder


    def arclength_to_r(self, s):
        t = self.arclength_to_theta(s)
        return self.theta_to_r(t)


    def get_arclength(self, vmin, vmax, var):
        '''
        Returns the arclength 
            for vmin <= theta <= vmax if var = 'theta'
            for vmin <= r <= vmax if var = 'r'

        '''
        if var not in ['r', 'theta']:
            raise ValueError('var (= "%s") undefined.'%var)
        if vmin < 0:
            raise ValueError('vmin (= %g) must be >= 0.'%vmin)
        if vmax < vmin:
            raise ValueError(
                'vmax (= %g) must be >= vmin (= %g).'%(vmax, vmin))
        if var == 'r':
            tmin = self.r_to_theta(vmin); tmax = self.r_to_theta(vmax)
        if var == 'theta':
            tmin = vmin; tmax = vmax

        tmp = tmax*math.sqrt(1.0+tmax*tmax) + math.asinh(tmax) - \
              tmin*math.sqrt(1.0+tmin*tmin) + math.asinh(tmin)
        return 0.5*self.b*tmp


    def get_curvature(self, v, var):
        '''
        Returns the curvature at a given value of theta, radius, or arclength.

        '''
        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        if var == 'r':
            t = self.r_to_theta(v)
        if var == 's':
            t = self.arclength_to_theta(v)
        if var == 'theta':
            t = v
        t_ = np.asarray(t)
        if np.any(t_ < self.tmin):
            raise ValueError('t must be >= lower bound (= %g).'%self.tmin)
        tsq = t_*t_
        den = self.b*(tsq +1)**1.5
        return (tsq+2)/den


    def get_avg_curv_rad(self, vmin, vmax, step, var):
        '''
        Returns the average radius of curvature over a given range of theta, radius, or arclength.

        '''
        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        if vmin < 0:
            raise ValueError('vmin (= "%g") must be >= 0.'%vmin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if vmin > vmax:
            raise ValueError(
                'vmin (= "%g") must be <= vmax (= "%g").'%(vmin, vmax))

        n = 1 + math.ceil( (vmax - vmin)/step )
        v = np.linspace(vmin, vmax, n)

        radii = np.zeros((n,))
        for i in range(n):
#           radii[i] = 1.0/self.get_curvature(v[i], var)
#       return radii.mean()
            radii[i] = self.get_curvature(v[i], var)
        return 1.0/radii.mean()


    def get_polar_coords(self, vmin, vmax, step, var, branch='out'):
        if branch not in ['in', 'out']:
            raise ValueError('branch (= "%s") undefined.'%branch)

        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        
        if vmin < 0:
            raise ValueError('vmin (= "%g") must be >= 0.'%vmin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if vmin > vmax:
            raise ValueError(
                'vmin (= "%g") must be <= vmax (= "%g").'%(vmin, vmax))

        n = 1 + math.ceil( (vmax - vmin)/step )
        v = np.linspace(vmin, vmax, n)

        if var == 'r':
            r = v; theta = self.r_to_theta(r)
        elif var == 'theta':
            theta = v; r = self.theta_to_r(theta)
        elif var == 's':
            theta = self.arclength_to_theta(v); r = self.theta_to_r(theta)

        if branch == 'in':
            theta += math.pi
        return r, theta


    def get_cart_coords(self, vmin, vmax, step, var, branch='out'):
        r, theta = self.get_polar_coords(vmin, vmax, step, var, branch)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return x, y




class DoubleSpiral(object):
    def __init__(self, spl, L, crvdir, is_tapered=False, taper_angle=0):
        if crvdir not in ['S', 'D']:
            raise ValueError('crvdir (= "%s") undefined.'%crvdir)
        self.spl = spl #A spiral object
        self.L = L
        self.crvdir = crvdir #Same or different curvature w.r.t. arclength
        self.is_tapered = is_tapered
        self.taper_angle = taper_angle

        if self.is_tapered:
            #Angle at the endpoint of one spiral
            te = self.spl.arclength_to_theta(self.L/2)
            #Angle at the beginning of taper
            t_tb = te - self.taper_angle
            #Length of the tapered arc
            self.taper_length = self.spl.get_arclength(t_tb, te, 'theta')
            #Arclength at the beginning of the tapered section
            self.spiral_length = self.L/2 - self.taper_length
            #Curvature at the beginning of the tapered section
            kappa = self.spl.get_curvature(self.spiral_length, 's')
            #Taper function
            tf = lambda x: kappa*(1 - x/self.taper_length)
            #Tapering spiral
            self.tspl = SpiralGeneral(tf)


    def get_curvature(self, s):
        '''
        Returns the curvature at a given value of arclength = s.

        '''
        if np.isscalar(s):
            s_ = np.array([s]); s_isscalar = True
        else:
            s_ = np.asarray(s); s_isscalar = False
        if np.any(s_ < 0):
            raise ValueError('s must be >= 0.')
        if np.any(s_ > self.L):
            raise ValueError('s must be <= total arclength (= %g).'%self.L)
        out = np.empty_like(s_)
        for i in range(s_.size):
            if s_[i] > self.L/2:
                s1 = self.L - s_[i]; ispr = 2
            else:
                s1 = s_[i]; ispr = 1

            if self.is_tapered and (s1 > self.spiral_length):
                kappa = self.tspl.get_curvature(s1-self.spiral_length)
            else:
                kappa = self.spl.get_curvature(s1, 's')

            if self.crvdir == 'D':
                out[i] = kappa if ispr == 1 else -kappa 
            elif self.crvdir == 'S':
                out[i] = kappa #if ispr == 1 else kappa 
        if s_isscalar:
            return out[0]
        else:
            return out


    def get_avg_curv_rad(self, vmin, vmax, step):
        '''
        Returns the average radius of curvature over a given range of arclength.

        '''
        if vmin < 0:
            raise ValueError('vmin (= "%g") must be >= 0.'%vmin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if vmin > vmax:
            raise ValueError(
                'vmin (= "%g") must be <= vmax (= "%g").'%(vmin, vmax))

        n = 1 + math.ceil( (vmax - vmin)/step )
        v = np.linspace(vmin, vmax, n)

        radii = np.zeros((n,))
        for i in range(n):
            radii[i] = self.get_curvature(v[i]) #Avoid infinite radius
        print(radii.mean())
        return 1.0/radii.mean()


    def get_cart_coords(self, step):
        if self.is_tapered:
            #Spiral section
            x1_s, y1_s = self.spl.get_cart_coords(0, self.spiral_length,
                            step, 's')
            #Tapering section
            x1_t, y1_t = self.tspl.get_cart_coords(step, self.taper_length, step)

            n = x1_s.size
            coords = np.empty((2,n))

            xe = x1_s[-1]; ye = y1_s[-1]
            ang = -(math.pi/2 + math.atan2(ye, xe))
            ct = math.cos(ang); st = math.sin(ang)
            rotmat = np.array([[ct, -st],[st, ct]])
            coords[0,:] = x1_s; coords[1,:] = y1_s
            tmp = np.dot(rotmat, coords)
            x1_s = tmp[0,:]
            y1_s = tmp[1,:] - tmp[1,-1]

            x1 = np.hstack((x1_s, x1_t))
            y1 = np.hstack((y1_s, y1_t))
        else:
            x1, y1 = self.spl.get_cart_coords(0, self.L/2, step, 's')
            n = x1.size
            coords = np.empty((2,n))
            xe = x1[-1]; ye = y1[-1]
            ang = -math.atan2(ye, xe)
            ct = math.cos(ang); st = math.sin(ang)
            rotmat = np.array([[ct, -st],[st, ct]])
            coords[0,:] = x1; coords[1,:] = y1
            tmp = np.dot(rotmat, coords)
            x1 = tmp[0,:]; y1 = tmp[1,:]

        n = x1.size
        coords = np.empty((2,2*n-1))
        coords[0,0:n] = x1; coords[1,0:n] = y1

        ce = coords[:,n-1]
        coords[0,0:n] -= ce[0]
        coords[1,0:n] -= ce[1]

        if self.crvdir == 'S':
            coords[0,n:] = -coords[0,n-1:0:-1]
            coords[1,n:] =  coords[1,n-1:0:-1]

        elif self.crvdir == 'D':
            #Flip for inner branch
            coords[0,n:] = -coords[0,n-1:0:-1]
            coords[1,n:] = -coords[1,n-1:0:-1]

        return coords[0,:], coords[1,:], n




class SpiralLogarithmic(object):
    '''
    Polar equation: r = a*exp(b*theta)

    '''
    def __init__(self, a, b, rmin):
        if a <= 0:
            raise ValueError('a (= %g) must be > 0.'%a)
        if b <= 0:
            raise ValueError('b (= %g) must be > 0.'%b)
        if rmin <= 0:
            raise ValueError('rmin (= %g) must be > 0.'%rmin)
        self.a = a; self.b = b; self.rmin = rmin
        self.tmin = self.r_to_theta(self.rmin)


    def r_to_theta(self, r):
        r_ = np.asarray(r)
        if np.any(r_ < self.rmin):
            raise ValueError('r must be >= rmin (= %g).'%self.rmin)
        return np.log(r_/self.a)/self.b


    def theta_to_r(self, theta):
        theta_ = np.asarray(theta)
        if np.any(theta_ < self.tmin):
            raise ValueError('theta must be >= tmin (= %g).'%self.tmin)
        return self.a*np.exp(self.b*theta_)


    def arclength_to_theta(self, s):
        '''
        Returns the theta parameter at arclength s. Note that here the 
        arclength s is measured from theta = tmin, where tmin is the starting
        value of theta for the spiral.

        '''
        r = self.arclength_to_r(s)
        return self.r_to_theta(r)


    def arclength_to_r(self, s):
        s_ = np.asarray(s)
        if np.any(s_ < 0):
            raise ValueError('s must be >= 0.')
        return self.rmin + self.b*s/math.sqrt(self.b*self.b+1)


    def get_arclength(self, vmin, vmax, var):
        '''
        Returns the arclength 
            for vmin <= theta <= vmax if var = 'theta'
            for vmin <= r <= vmax if var = 'r'

        '''
        if var not in ['r', 'theta']:
            raise ValueError('var (= "%s") undefined.'%var)
        if vmin < 0:
            raise ValueError('vmin (= %g) must be >= 0.'%vmin)
        if vmax < vmin:
            raise ValueError(
                'vmax (= %g) must be >= vmin (= %g).'%(vmax, vmin))
        if var == 'r':
            tmp = vmax - vmin
        elif var == 'theta':
            tmp = self.theta_to_r(vmax) - self.theta_to_r(vmin)
        return math.sqrt(self.b*self.b+1)*tmp/self.b


    def get_curvature(self, v, var):
        '''
        Returns the curvature at a given value of theta, radius, or arclength.

        '''
        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        if var == 'r':
            r = v
        if var == 's':
            r = self.arclength_to_r(v)
        if var == 'theta':
            r = self.theta_to_r(v)
        r_ = np.asarray(r)
        if np.any(r_ < self.rmin):
            raise ValueError('r must be >= lower bound (= %g).'%self.rmin)
        tmp = r_*math.sqrt(1+self.b*self.b)
        return 1.0/tmp


    def get_polar_coords(self, vmin, vmax, step, var):
        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        
        if vmin < 0:
            raise ValueError('vmin (= "%g") must be >= 0.'%vmin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if vmin > vmax:
            raise ValueError(
                'vmin (= "%g") must be <= vmax (= "%g").'%(vmin, vmax))

        n = 1 + math.ceil( (vmax - vmin)/step )
        v = np.linspace(vmin, vmax, n)

        if var == 'r':
            r = v; theta = self.r_to_theta(r)
        elif var == 'theta':
            theta = v; r = self.theta_to_r(theta)
        elif var == 's':
            r = self.arclength_to_r(v);
            theta = self.r_to_theta(r)
        return r, theta


    def get_cart_coords(self, vmin, vmax, step, var):
        r, theta = self.get_polar_coords(vmin, vmax, step, var)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return x, y




class SpiralFermat(object):
    '''
    Polar equation: r = b*theta^(1/2)

    '''
    def __init__(self, b, rmin):
        if b <= 0:
            raise ValueError('b (= %g) must be > 0.'%b)
        if rmin < 0:
            raise ValueError('rmin (= %g) must be >= 0.'%rmin)
        self.b = b
        self.rmin = rmin
        self.tmin = self.r_to_theta(self.rmin)


    def r_to_theta(self, r):
        r_ = np.asarray(r)
        if np.any(r_ < self.rmin):
            raise ValueError('r must be >= rmin (= %g).'%self.rmin)
        return (r_/self.b)**2


    def theta_to_r(self, theta):
        theta_ = np.asarray(theta)
        if np.any(theta_ < self.tmin):
            raise ValueError('theta must be >= tmin (= %g).'%self.tmin)
        return self.b*np.sqrt(theta_)


    def arclength_to_r(self, s):
        '''
        Returns the radius at arclength s. Note that here the 
        arclength s is measured from r = rmin, where rmin is the starting
        value of r for the spiral.

        '''
        if np.isscalar(s):
            s_ = np.array([s]); s_isscalar = True
        else:
            s_ = np.asarray(s); s_isscalar = False
        if np.any(s_ < 0):
            raise ValueError('s must be >= 0.')
        out = np.empty_like(s_)
        n = s_.size
        for i in range(n):
            sol = root_scalar(self.__class__._al2r_func, 
                    args=(self.b, self.rmin, s_[i]), 
                    method='newton', fprime=True, x0=self.rmin, xtol=1.e-8,
                    rtol=1.e-10)
            #Check convergence
            if not sol.converged:
                print('r calculation not converged for element %d'%i)
            out[i] = sol.root
        if s_isscalar:
            return out[0]
        else:
            return out


    @staticmethod
    def _al2r_func(r, *args):
        b = args[0]; rmin = args[1]; s = args[2] 
        rmin4 = (rmin/b)**4; r4 = (r/b)**4
        C = rmin*hyp2f1(-0.5, 0.25, 1.25, -4*rmin4)
        f = r*hyp2f1(-0.5, 0.25, 1.25, -4*r4) - (s + C)
        #See https://dlmf.nist.gov/15.5  (Eq 15.5.1)
        fder = hyp2f1(-0.5, 0.25, 1.25, -4*r4) + \
               2*r4*hyp2f1(0.5, 1.25, 2.25, -4*r4)
        return f, fder


    def arclength_to_theta(self, s):
        r = self.arclength_to_r(s)
        return self.r_to_theta(r)


    def get_arclength(self, vmin, vmax, var):
        '''
        Returns the arclength 
            for vmin <= theta <= vmax if var = 'theta'
            for vmin <= r <= vmax if var = 'r'

        '''
        if var not in ['r', 'theta']:
            raise ValueError('var (= "%s") undefined.'%var)
        if vmin < 0:
            raise ValueError('vmin (= %g) must be >= 0.'%vmin)
        if vmax < vmin:
            raise ValueError(
                'vmax (= %g) must be >= vmin (= %g).'%(vmax, vmin))
        if var == 'r':
            tmin = self.r_to_theta(vmin); tmax = self.r_to_theta(vmax)
        if var == 'theta':
            tmin = vmin; tmax = vmax

        tmp = math.sqrt(tmax)*hyp2f1(-0.5, 0.25, 1.25, -4*tmax*tmax) - \
              math.sqrt(tmin)*hyp2f1(-0.5, 0.25, 1.25, -4*tmin*tmin)
        return self.b*tmp


    def get_curvature(self, v, var):
        '''
        Returns the curvature at a given value of theta, radius, or arclength.

        '''
        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        if var == 'r':
            r = v
        if var == 's':
            r = self.arclength_to_r(v)
        if var == 'theta':
            r = self.theta_to_r(v)
        r_ = np.asarray(r)
        if np.any(r_ < self.rmin):
            raise ValueError('r must be >= lower bound (= %g).'%self.rmin)
        r4 = r_**4
        return 2*r_*(4*r4 + 3*self.b**4)/(4*r4 + self.b**2)**1.5


    def get_polar_coords(self, vmin, vmax, step, var, branch='out'):
        if branch not in ['in', 'out']:
            raise ValueError('branch (= "%s") undefined.'%branch)

        if var not in ['r', 'theta', 's']:
            raise ValueError('var (= "%s") undefined.'%var)
        
        if vmin < 0:
            raise ValueError('vmin (= "%g") must be >= 0.'%vmin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if vmin > vmax:
            raise ValueError(
                'vmin (= "%g") must be <= vmax (= "%g").'%(vmin, vmax))

        n = 1 + math.ceil( (vmax - vmin)/step )
        v = np.linspace(vmin, vmax, n)

        if var == 'r':
            r = v; theta = self.r_to_theta(r)
        elif var == 'theta':
            theta = v; r = self.theta_to_r(theta)
        elif var == 's':
            theta = self.arclength_to_theta(v); r = self.theta_to_r(theta)

        if branch == 'in':
            theta += math.pi
        return r, theta


    def get_cart_coords(self, vmin, vmax, step, var, branch='out'):
        r, theta = self.get_polar_coords(vmin, vmax, step, var, branch)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return x, y




class SpiralClothoid(object):
    '''
    Ref. : https://mathcurve.com/courbes2d.gb/cornu/cornu.shtml
           https://mathworld.wolfram.com/CornuSpiral.html

    Cartesian equation: x = a*C(t); y = a*S(t), where C(t) and S(t) are Frenel
    integrals. Note that the argument of the Frenel integerals is (pi/2)*t^2
    rather than just t^2. Asymptotic point is (a/2, a/2).

    '''
    def __init__(self, a):
        if a <= 0:
            raise ValueError('a (= %g) must be > 0.'%a)
        self.a = a


    def get_curvature(s): 
        return s/(self.a*self.a)


    def get_cart_coords(self, smin, smax, step):
        '''
        Arclength is measured from the point of zero curvature.
        '''
        if smin < 0:
            raise ValueError('smin (= "%g") must be >= 0.'%smin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if smin > smax:
            raise ValueError(
                'smin (= "%g") must be <= smax (= "%g").'%(smin, smax))

        n = 1 + math.ceil( (smax - smin)/step )
        s = np.linspace(smin/self.a, smax/self.a, n)
        y, x = fresnel(s)
        x *= self.a; y *= self.a
        return x, y




class SpiralGeneral(object):
    '''
    Cartesian equation: x = int_0^s cos(int_kappa) ds
                        y = int_0^s sin(int_kappa) ds
    where int_kappa = int_0^t kappa dt

    '''
    def __init__(self, func_curvature):
        self.func_kappa = func_curvature


    def _int_kappa(self, x):
        out = np.empty_like(x)
        for i in range(out.shape[0]):
            out[i], err = quadrature(self.func_kappa, 0, x[i], maxiter=200)
        return out


    def get_curvature(self, s):
        return self.func_kappa(s)


    def get_cart_coords(self, smin, smax, step):
        if smin < 0:
            raise ValueError('smin (= "%g") must be >= 0.'%smin)
        if step <= 0:
            raise ValueError('step (= "%g") must be > 0.'%step)
        if smin > smax:
            raise ValueError(
                'smin (= "%g") must be <= smax (= "%g").'%(smin, smax))

        n = 1 + math.ceil( (smax - smin)/step )
        s = np.linspace(smin, smax, n)
        x = np.empty_like(s); y = np.empty_like(s)

        f = lambda x: np.cos(self._int_kappa(x))
        g = lambda x: np.sin(self._int_kappa(x))
        for i in range(n):
            x[i],err = quadrature(f, 0, s[i], maxiter=400)
            y[i],err = quadrature(g, 0, s[i], maxiter=400)
        return x, y


    def get_polar_coords(self, smin, smax, step):
        x, y = self.get_cart_coords(smin, smax, step)
        r = np.hypot(x, y)
        theta = np.atan2(y, x)
        return r, theta






