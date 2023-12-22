#!/usr/bin/env python

'''
Generates a monolayer sheet with a spiral profile

'''

import math
import numpy as np
from spirals import *

def write_xyz(fn, fn_title, xdata, ydata, zdata, nx, nz):
    nnodes = nx*nz
    nodes = np.zeros((nnodes,3))
    for k in range(nz):
        for i in range(nx):
            inode = k*nx + i
            nodes[inode,0] = xdata[i]
            nodes[inode,1] = ydata[i]
            nodes[inode,2] = zdata[k]
    
    with open(fn, 'w') as fh:
        fh.write(str(nnodes) + '\n')
        fh.write(fn_title + '\n')
        for i in range(nnodes):
            fh.write( 'M  ' + '  '.join(['% .15E'%x for x in nodes[i,:]]) + '\n')


#Sheet dimensions
Lx = 100 #Along x
Lz = 10  #Along z

dx = 0.5
dz = 0.5

nx = int(Lx/dx)
nz = int(Lz/dz)

zdata = np.linspace(0, Lz, nz)

#Generate profile curves

#Spiral archimedean
#
cnst_b = 0.47 # 0.2, 0.3, 0.4, 0.47, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6
rmin = 3.62 #3.62 [0, 0.8, 1.6, 2.4, 3, 3.62, 4, 5]

for cnst_b in [0.47]:
#   for rmin in [0, 0.8, 1.6, 2.4, 3.0, 3.62, 4.0, 5.0]:
        spl = SpiralArchimedes(cnst_b, rmin)
        xdata, ydata = spl.get_cart_coords(0, Lx, dx, var='s')
        fn = 'sa_b_%g_rmin_%g_L_%g.xyz'%(cnst_b, rmin, Lx)
        fn_title = 'spiral_arch'
        write_xyz(fn, fn_title, xdata, ydata, zdata, nx, nz)
        
        #taper_angle = 35 for b >= 1.2
        #taper_angle = 50 for 0.2 <= b < 1.2
        #ta = 45
        #ds = DoubleSpiral(spl, Lx, 'D', is_tapered=True, taper_angle=math.radians(ta))
        #xdata, ydata, tmp = ds.get_cart_coords(dx)
        #fn = 'dsa_D_b_%g_rmin_%g_L_%g_ta_%g.xyz'%(cnst_b, rmin, Lx, ta)
        #fn_title = 'double_spiral_arch'
        #write_xyz(fn, fn_title, xdata, ydata, zdata, nx, nz)

#Spiral fermat
#
#cnst_b = 1.6; rmin = 0.0
 
#spl = SpiralFermat(cnst_b, rmin)
#xdata, ydata = spl.get_cart_coords(0, Lx, dx, var='s')
#fn = 'sf_b_%g_rmin_%g.xyz'%(cnst_b, rmin)
#fn_title = 'spiral_fermat'
 
#ds = DoubleSpiral(spl, Lx, 'S', is_tapered=True, taper_angle=math.radians(35))
#xdata, ydata = ds.get_cart_coords(0, Lx, dx)
#fn = 'dsf_D_b_%g_rmin_%g.xyz'%(cnst_b, rmin)
#fn_title = 'double_spiral_fermat'

#Spiral logarithmic
#cnst_a = 1.0; cnst_b = 0.2; rmin = 1.0

#spl = SpiralLogarithmic(cnst_a, cnst_b, rmin)
#xdata, ydata = spl.get_cart_coords(0, Lx, dx, var='s')
#fn = 'sl_a_%g_b_%g_rmin_%g.xyz'%(cnst_a, cnst_b, rmin)
#fn_title = 'spiral_logarithmic'
 
#ds = DoubleSpiral(spl, Lx, 'S', is_tapered=True, taper_angle=math.radians(35))
#xdata, ydata = ds.get_cart_coords(0, Lx, dx)
#fn = 'dsl_S_a_%g_b_%g_rmin_%g.xyz'%(cnst_a, cnst_b, rmin)
#fn_title = 'double_spiral_logarithmic'


#Spiral clothoid
#cnst_a = 400.0
 
#spl = SpiralClothoid(cnst_a)
#xdata, ydata = spl.get_cart_coords(0, Lx, dx)
#fn = 'sc_a_%g.xyz'%(cnst_a)
#fn_title = 'spiral_clothoid'
 
#x1, y1 = spl.get_cart_coords(0, Lx, dx)
#x2 = -x1; y2 = y1   # Curvdir = 'S'
#x2 = -x1; y2 = -y1  # Curvdir = 'D'
#xdata = np.hstack((x1, x2))
#ydata = np.hstack((y1, y2))
#fn = 'dsc_S_a_%g_b_%g_rmin_%g.xyz'%(cnst_a, rmin)
#fn_title = 'double_spiral_clothoid'
    
    
#Plane
#xdata = np.linspace(0, Lx, nx)
#ydata = np.zeros((nx,))
#fn = 'flat_%dx%d.xyz'%(Lx, Lz)
#fn_title = 'plane %dx%d'%(Lx, Lz)
#write_xyz(fn, fn_title, xdata, ydata, zdata, nx, nz)
