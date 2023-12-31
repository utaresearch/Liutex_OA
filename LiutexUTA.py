#!/bin/python3
##########################################################################
###  Liutex module for useful functions in python.
###  
###  Author: Oscar Alvarez
###  Email_1: oscar.alvarez@uta.edu
###  Email_2: mathOscarAlvarez@uta.edu
###  Center for Numerical Simulation and Modeling
###  University of Texas at Arlington
###  Department of Mathematics
##########################################################################
import numpy as np

def velocity_gradient3d(x : list, y : list, z : list, u : list, v : list, w : list) -> list:
    ## Calculates the velocity gradient tensor using the u, v, and w velocities using
    ## the global coordinates x, y, and z.

    imax = x.shape[0]
    jmax = x.shape[1]
    kmax = x.shape[2]

    velocity_gradient = np.empty([imax, jmax, kmax, 3, 3])

    for k in range(kmax):
        for j in range(jmax):
            for i in range(imax):
                ## Using Finite difference scheme to calculate partial derivatives

                if (i == 0):
                    ## forward difference
                    u_xi = u[1, j, k] - u[0, j, k]
                    v_xi = v[1, j, k] - v[0, j, k]
                    w_xi = w[1, j, k] - w[0, j, k]
        
                    x_xi = x[1, j, k] - x[0, j, k]
                    y_xi = y[1, j, k] - y[0, j, k]
                    z_xi = z[1, j, k] - z[0, j, k]
                elif (i == imax-1):
                    ## backward difference
                    u_xi = u[imax-1, j, k] - u[imax-2, j, k]
                    v_xi = v[imax-1, j, k] - v[imax-2, j, k]
                    w_xi = w[imax-1, j, k] - w[imax-2, j, k]
        
                    x_xi = x[imax-1, j, k] - x[imax-2, j, k]
                    y_xi = y[imax-1, j, k] - y[imax-2, j, k]
                    z_xi = z[imax-1, j, k] - z[imax-2, j, k]
                else:
                    ## central difference
                    u_xi = 0.5*(u[i+1, j, k] - u[i-1, j, k])
                    v_xi = 0.5*(v[i+1, j, k] - v[i-1, j, k])
                    w_xi = 0.5*(w[i+1, j, k] - w[i-1, j, k])
        
                    x_xi = 0.5*(x[i+1, j, k] - x[i-1, j, k])
                    y_xi = 0.5*(y[i+1, j, k] - y[i-1, j, k])
                    z_xi = 0.5*(z[i+1, j, k] - z[i-1, j, k])
    

                if (j == 0):
                    ## forward difference
                    u_eta = u[i, 1, k] - u[i, 0, k]
                    v_eta = v[i, 1, k] - v[i, 0, k]
                    w_eta = w[i, 1, k] - w[i, 0, k]
        
                    x_eta = x[i, 1, k] - x[i, 0, k]
                    y_eta = y[i, 1, k] - y[i, 0, k]
                    z_eta = z[i, 1, k] - z[i, 0, k]
                elif (j == jmax-1):
                    ## backward difference
                    u_eta = u[i, jmax-1, k] - u[i, jmax-2, k]
                    v_eta = v[i, jmax-1, k] - v[i, jmax-2, k]
                    w_eta = w[i, jmax-1, k] - w[i, jmax-2, k]
        
                    x_eta = x[i, jmax-1, k] - x[i, jmax-2, k]
                    y_eta = y[i, jmax-1, k] - y[i, jmax-2, k]
                    z_eta = z[i, jmax-1, k] - z[i, jmax-2, k]    
                else:
                    ## central difference
                    u_eta = 0.5*(u[i, j+1, k] - u[i, j-1, k])
                    v_eta = 0.5*(v[i, j+1, k] - v[i, j-1, k])
                    w_eta = 0.5*(w[i, j+1, k] - w[i, j-1, k])
        
                    x_eta = 0.5*(x[i, j+1, k] - x[i, j-1, k])
                    y_eta = 0.5*(y[i, j+1, k] - y[i, j-1, k])
                    z_eta = 0.5*(z[i, j+1, k] - z[i, j-1, k])    
                
    
                if (k == 0):    
                    ## forward difference
                    u_zeta = u[i, j, 1] - u[i, j, 0]
                    v_zeta = v[i, j, 1] - v[i, j, 0]
                    w_zeta = w[i, j, 1] - w[i, j, 0]
        
                    x_zeta = x[i, j, 1] - x[i, j, 0]
                    y_zeta = y[i, j, 1] - y[i, j, 0]
                    z_zeta = z[i, j, 1] - z[i, j, 0]
                elif (k == kmax-1):
                    ## backward difference
                    u_zeta = u[i, j, kmax-1] - u[i, j, kmax-2]
                    v_zeta = v[i, j, kmax-1] - v[i, j, kmax-2]
                    w_zeta = w[i, j, kmax-1] - w[i, j, kmax-2]
        
                    x_zeta = x[i, j, kmax-1] - x[i, j, kmax-2]
                    y_zeta = y[i, j, kmax-1] - y[i, j, kmax-2]
                    z_zeta = z[i, j, kmax-1] - z[i, j, kmax-2]
                else:
                    ## central difference
                    u_zeta = 0.5*(u[i, j, k+1] - u[i, j, k-1])
                    v_zeta = 0.5*(v[i, j, k+1] - v[i, j, k-1])
                    w_zeta = 0.5*(w[i, j, k+1] - w[i, j, k-1])
        
                    x_zeta = 0.5*(x[i, j, k+1] - x[i, j, k-1])
                    y_zeta = 0.5*(y[i, j, k+1] - y[i, j, k-1])
                    z_zeta = 0.5*(z[i, j, k+1] - z[i, j, k-1])

    
                ## determinant of Jacobian
                det = x_xi * (y_eta*z_zeta-y_zeta*z_eta) - x_eta * (y_xi*z_zeta-y_zeta*z_xi) + x_zeta * (y_xi*z_eta-y_eta*z_xi)
                
                if (det == 0):
                    det = 0.0
                else:
                    det = 1.0 / det
    
                xi_x = det*(y_eta*z_zeta - y_zeta*z_eta)
                xi_y = det*(x_zeta*z_eta - x_eta*z_zeta)
                xi_z = det*(x_eta*y_zeta - x_zeta*y_eta)
    
                eta_x = det*(y_zeta*z_xi - y_xi*z_zeta)
                eta_y = det*(x_xi*z_zeta - x_zeta*z_xi)
                eta_z = det*(x_zeta*y_xi - x_xi*y_zeta)
    
                zeta_x = det*(y_xi*z_eta - y_eta*z_xi)
                zeta_y = det*(x_eta*z_xi - x_xi*z_eta)
                zeta_z = det*(x_xi*y_eta - x_eta*y_xi)

                ## Assembling the velocity gradient tensor
                
                velocity_gradient[i, j, k, 0, 0] = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x    ##dudx
                velocity_gradient[i, j, k, 0, 1] = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y    ##dudy
                velocity_gradient[i, j, k, 0, 2] = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z    ##dudz

                velocity_gradient[i, j, k, 1, 0] = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x    ##dvdx
                velocity_gradient[i, j, k, 1, 1] = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y    ##dvdy
                velocity_gradient[i, j, k, 1, 2] = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z    ##dvdz
                
                velocity_gradient[i, j, k, 2, 0] = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x    ##dwdx
                velocity_gradient[i, j, k, 2, 1] = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y    ##dwdy
                velocity_gradient[i, j, k, 2, 2] = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z    ##dwdz
            
    return velocity_gradient


def liutex_3d(velocity_gradient_matrix : list) -> tuple([float, list]):
    """Calculates the 3D Liutex magnitude and Liutex direction vector using
    the velocity gradient tensor (matrix).
    Returns: R, r (Liutex Magnitude, Liutex vector)
    University of Texas at Arlington - Department of Mathematics
    Author: Oscar Alvarez
    Email: oscar.alvarez@uta.edu
    """
    # velocity gradient matrix = [[du/dx du/dy du/dz],
    #                             [dv/dx dv/dy dv/dz],
    #                             [dw/dx dw/dy dw/dz]]

    A = np.array(velocity_gradient_matrix)

    eig_val, eig_vec = np.linalg.eig(A)

    # If no complex eigenvalues exist, there is no rotation and Liutex = 0
    if np.isrealobj(eig_val):
        return 0, np.array([0.0, 0.0, 0.0])
    
    # Find location (index) of real eigenvalue/eigenvector
    for i in range(3):
        try:
            if (eig_val[i].imag == 0.0):
                real_index = i
                complex_index = (i+1) % 3
                break
        except:
            print('i = ', i)
            print('eig_val shape: ', eig_val[i].shape)
            exit()


    real_eig_vec = eig_vec[:,real_index].real
    
    # Complex Eigenvalue
    lambda_ci = eig_val[complex_index].imag
    
    # Find vorticity
    vorticity = np.zeros(3)
    vorticity[0] = A[2,1] - A[1,2]
    vorticity[1] = A[0,2] - A[2,0]
    vorticity[2] = A[1,0] - A[0,1]

    # Liutex direction condition
    if (np.dot(vorticity, real_eig_vec) < 0):
        r = -real_eig_vec.copy()
    else:
        r = real_eig_vec.copy()

    # Finding Liutex Magnitude
    w_dot_r = np.dot(vorticity, r)
    liutex_magnitude = w_dot_r - np.sqrt( w_dot_r*w_dot_r - 4.0*lambda_ci*lambda_ci )

    # Liutex vector r
    r = liutex_magnitude * r

    return liutex_magnitude, r


def modified_omega_liutex_3d(velocity_gradient_matrix : list) -> tuple([float, list]):
    """Calculates the 3D Modified Omega Liutex magnitude and vector using
    the velocity gradient tensor (matrix).
    University of Texas at Arlington - Department of Mathematics
    Author: Oscar Alvarez
    Email: oscar.alvarez@uta.edu
    """
    # velocity gradient matrix = [[du/dx du/dy du/dz],
    #                             [dv/dx dv/dy dv/dz],
    #                             [dw/dx dw/dy dw/dz]]

    A = np.array(velocity_gradient_matrix)

    eig_val, eig_vec = np.linalg.eig(A)

    # If no complex eigenvalues exist, there is no rotation and Liutex = 0
    if np.isrealobj(eig_val):
        return 0, np.array([0.0, 0.0, 0.0])
    
    # Find real eigenvalue/eigenvector
    for i in range(3):
        if eig_val[i].imag == 0.0:
            real_index = np.copy(i)
            complex_index = (i+1) % 3
            break
    
    lambda_r = eig_val[real_index].real
    real_eig_vec = eig_vec[:,real_index].real
    
    # Complex Eigenvalue
    lambda_cr = eig_val[complex_index].real
    lambda_ci = eig_val[complex_index].imag
    
    # Find vorticity
    vorticity = np.zeros(3)
    vorticity[0] = A[2,1] - A[1,2]
    vorticity[1] = A[0,2] - A[2,0]
    vorticity[2] = A[1,0] - A[0,1]

    # Liutex direction condition
    if (np.dot(vorticity, real_eig_vec) < 0):
        r = -real_eig_vec.copy()
    else:
        r = real_eig_vec.copy()

    epsilon = 1.0e-3

    # Finding Modified Omega Liutex Magnitude
    w_dot_r = np.dot(vorticity, r)
    modified_omega_liutex = w_dot_r*w_dot_r / (2.0*(w_dot_r*w_dot_r - 2.0*lambda_ci*lambda_ci + 2.0*lambda_cr*lambda_cr + lambda_r*lambda_r) + epsilon)

    # Modified Omega Liutex vector
    Or = modified_omega_liutex * r

    return modified_omega_liutex, Or


def write_PLOT3D(filename_without_extension : str, grid_data: list, function_data : list):
    """Creates PLOT3D files (.xyz) and (.fun) from grid and function data.
    Grid data format:
    grid_data = np.empty([imax, jmax, kmax, 3])
    grid_data[:,:,:,0] = x[:,:,:]
    grid_data[:,:,:,1] = y[:,:,:]
    grid_data[:,:,:,2] = z[:,:,:]
    """

    print('file: ', filename_without_extension)

    imax = grid_data.shape[0]
    jmax = grid_data.shape[1]
    kmax = grid_data.shape[2]

    grid_dim_list = [str(imax), str(jmax), str(kmax)]
    grid_dim_string = ' '.join(grid_dim_list)

    out_grid_filename = filename_without_extension + '_grid.xyz'

    print('Writing Grid File (.xyz)')
    
    with open(out_grid_filename, 'w') as grid_f:
        grid_f.write(grid_dim_string + '\n')

        for n in range(3):
            for k in range(kmax):
                for j in range(jmax):
                    for i in range(imax):
                        grid_f.write(str(grid_data[i,j,k,n]) + ' ')


    print('Writing Function File (.fun)')
    
    nvars = function_data.shape[3]

    function_dim_list = [str(imax), str(jmax), str(kmax), str(nvars)]
    function_dim_string = ' '.join(function_dim_list)

    out_function_filename = filename_without_extension + '_function.fun'

    with open(out_function_filename, 'w') as function_f:
        function_f.write(function_dim_string + '\n')

        for n in range(nvars):
            for k in range(kmax):
                for j in range(jmax):
                    for i in range(imax):
                        function_f.write(str(function_data[i,j,k,n]) + ' ')
    
    print('PLOT3D files written successfully.')


if __name__ == "__main__":
    """Example use of liutex_3d() using example velocity gradient matrix.
    University of Texas at Arlington - Department of Mathematics
    Author: Oscar Alvarez
    Email: oscar.alvarez@uta.edu
    """
    example_vel_grad_matrix = [[4.0, -3.0,  0.0],
                               [3.0,  4.0,  0.0],
                               [5.0, 10.0, 10.0]]
    
    R, r = liutex_3d(example_vel_grad_matrix)

    print("\nLiutex Magnitude (R): ", R)
    print('Liutex Vector (r): ', r)

    OR, Or = modified_omega_liutex_3d(example_vel_grad_matrix)

    print("\nModified Omega Liutex Magnitude (OR): ", OR)
    print('Modified Omega Liutex Liutex Vector (Or): ', Or)
    