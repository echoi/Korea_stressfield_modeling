import numpy as np
from scipy import linalg as la

s_xx = np.array([-150.0])
s_yy = np.array([-50.0])
s_zz = np.array([-100.0])
s_xy = np.array([0.0])
s_yz = np.array([0.0])
s_xz = np.array([0.0])

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
  print(u,indx,sigma_d)
  print(v[0,:])
  print(v[1,:])
  print(v[2,:])

  sigma_3.append(v[indx[2]].tolist())
  sigma_1.append(v[indx[0]].tolist())
  sigma_2.append(v[indx[1]].tolist())

  mag_3.append(u[indx[2]].tolist())
  mag_1.append(u[indx[0]].tolist())
  mag_2.append(u[indx[1]].tolist())

# strike slip component 10 degree left

dips = np.linspace(10.0, 100.0, 1)
strikes = np.linspace(0.0, 370.0, 1)
for k in range(len(dips)):
  for j in range(len(strikes)):
    theta = np.radians(strikes[j]) # strike in radian
    coloumb = []
    phi = np.radians(40) ; # along z
    gamma = np.radians(dips[k]); # along x # dip in radian
    print(k,j,'strike=',theta,' dip=',gamma)
    # rotation matrix for rotation about z
    R = np.array([[np.cos(theta), -np.sin(theta), 0 ], [np.sin(theta), np.cos(theta), 0 ], [0, 0, 1]])
    Rx_dip = np.array([[1, 0, 0], [0, np.cos(gamma), -np.sin(gamma)], [0, np.sin(gamma), np.cos(gamma) ]])

    Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi) ]])
    
    print(R, Rx_dip, Rx)
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
    #output.PointData.append(np.array(coloumb), "coloumb_s%s_d%s" %(j, k))