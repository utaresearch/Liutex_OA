import numpy as np

# Data grid dimensions
ni = 20
nj = 20
nk = 50    
step_size = 0.5      # Step size

# Function constants
xi = 0.042
gma = 1.45
re = 23.1
nu = gma/(2.0*np.pi*re)


with open('burger_data.dat', 'w') as f:

    f.write(r'variables ="x","y","z","u","v","w","du_dx","du_dy","du_dz","dv_dx","dv_dy","dv_dz","dw_dx","dw_dy","dw_dz"')
    
    print('zone f=point, I=', ni+1, ', J=', str(nj+1), ', K=', str(nk+1), sep='   ', file=f )

    for k in range(nk+1):
        for j in range(nj+1):
            for i in range(ni+1):

                x = step_size*(i - ni/2.0)
                y = step_size*(j - nj/2.0)
                z = step_size*(k - nk/2.0)
				
                if (x == 0.0 and y == 0.0):
                    radi = 1e-7            # just to avoid NAN
                else:
                    radi = np.sqrt(x**2 + y**2)
				
				# velocity components of burger data (be careful while choosing parameter;vortex fully depends on parameters.)
                u = -xi*x - (gma/(2.0*np.pi*radi*radi))*(1.0 - np.exp(-radi*radi*xi/(2.0*nu)))*y
                v = -xi*y + (gma/(2.0*np.pi*radi*radi))*(1.0 - np.exp(-radi*radi*xi/(2.0*nu)))*x
                w = 2.0*xi*z
								
                # Partial Derivatives
                u_x = -xi
                u_y = -(gma/(2.0*np.pi*radi*radi))*(1.0 - np.exp(-radi*radi*xi/(2.0*nu)))
                u_z = 0.0
				
                v_x = (gma/(2.0*np.pi*radi*radi))*(1.0 - np.exp(-radi*radi*xi/(2.0*nu)))
                v_y = -xi
                v_z = 0.0
				
                w_x = 0.0
                w_y = 0.0
                w_z = 2.0*xi

                # Write values to file
                print(x, y, z, u, v, w, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, sep='   ', file=f)
