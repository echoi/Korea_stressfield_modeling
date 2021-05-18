#!/usr/bin/env pvpython
import numpy as np
from scipy import linalg
from vtk.numpy_interface import dataset_adapter as da

s_xx = inputs[0].PointData['stress_xx']
s_yy = inputs[0].PointData['stress_yy']
s_zz = inputs[0].PointData['stress_zz']
s_xy = inputs[0].PointData['stress_xy']
s_yz = inputs[0].PointData['stress_yz']
s_xz = inputs[0].PointData['stress_xz']

direction = []
sigma_1 = []
sigma_d = []; mag_1 = [];
sigma_2 = []; mag_2 = [];
sigma_3 = []; mag_3 = [];

for i in range(len(s_xx)):
  stress = np.array ( [[s_xx[i], s_xy[i], \
  s_xz[i]], [s_xy[i], s_yy[i], s_yz[i] ], [s_xz[i], \
  s_yz[i], s_zz[i] ]] )
  u,v  = np.linalg.eig(stress)
  indx = np.argsort(u)
  sigma_d.append(u[indx[2]] - u[indx[0]]) 

  sigma_3.append(v[indx[2]].tolist())
  sigma_1.append(v[indx[0]].tolist())
  sigma_2.append(v[indx[1]].tolist())

  mag_3.append(u[indx[2]].tolist())
  mag_1.append(u[indx[0]].tolist())
  mag_2.append(u[indx[1]].tolist())

vtk_arr = da.VTKArray(sigma_3)
vtk_arr2 = da.VTKArray(sigma_2)
vtk_arr3 = da.VTKArray(sigma_1)

output.PointData.append(vtk_arr, "sigma_1")
output.PointData.append(vtk_arr2, "sigma_2")
output.PointData.append(vtk_arr3, "sigma_3")

output.PointData.append(np.array(mag_3), "mag_3")
output.PointData.append(np.array(mag_2), "mag_2")
output.PointData.append(np.array(mag_1), "mag_1")

output.PointData.append(np.array(sigma_d), "sigma_d")

####### routine for stress rotation and coloumb stress calculation

# strike slip component 10 degree left
for k in range(10, 100, 10):
  for j in range(0, 370, 10):
    theta = np.radians(j)
    coloumb = []
    phi = np.radians(40) ; # along z
    gamma = np.radians(k); # along x

    # rotation matrix for rotation about z
    R = [[np.cos(theta), -np.sin(theta), 0 ], [np.sin(theta), np.cos(theta), 0 ], [0, 0, 1]]
    Rx_dip = [[1, 0, 0], [0, np.cos(gamma), -np.sin(gamma)], [0, np.sin(gamma), np.cos(gamma) ]]

    Rx = [[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi) ]]
    
    for i in range(len(s_xx)):
      stress = np.array ( [[s_xx[i], s_xy[i], s_xz[i]], [ s_xy[i], \
        s_yy[i], s_yz[i] ], [s_xz[i], s_yz[i], s_zz[i] ] ] )
      
      stress_dash_x = np.dot (np.dot (np.transpose(Rx), stress), Rx)
      stress_dash = np.dot (np.dot (np.transpose(R), stress_dash_x), R) # fault plane
      stress_dash_xx = np.dot (np.dot (np.transpose(Rx_dip), stress_dash), Rx_dip)
      
      sigma_nn = stress_dash_xx[2][2] # normal stress sigma_yy 
      tau = stress_dash_xx[2][1] # assuming strike-slip component sigma_xy
      mu  = 0.6 # coefficient of friction
      cff = tau + mu*sigma_nn
      coloumb.append(cff)
    output.PointData.append(np.array(coloumb), "coloumb_s%s_d%s" %(j, k))