import numpy as np


def liutex_3d(velocity_gradient_matrix : list) -> tuple([float, list]):
    """Calculates the 3D Liutex magnitude and Liutex direction vector using
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
    
    # Find location (index) of real eigenvalue/eigenvector
    for i in range(3):
        if eig_val[i].imag == 0.0:
            real_index = i
            complex_index = (i+1) % 3
            break

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
    