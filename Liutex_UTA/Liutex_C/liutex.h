////////////////////////////////////////////////////////////////////////////////////
/// Liutex method C code.
/// This code calculates liutex when given the veloctiy gradient tensor 
/// for 3D flow fields.
/// 
/// Author:  Oscar Alvarez
/// 
/// email: oscar.alvarez@uta.edu
/// University of Texas at Arlington
/// Department of Mathematics
/// Center for Numerical Simulation and Modeling (CNSM)
/// Arlington, Texas, United States
////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>



// Index functions for data read and manipulation.
int get_index3(int i, int j, int k, 
               int imax, int jmax, int kmax)
{
    return i + j*imax + k*imax*jmax;
}


int get_index4(int i, int j, int k, int n, 
               int imax, int jmax, int kmax, int nvar)
{
    return i + j*imax + k*imax*jmax + n*imax*jmax*kmax;
}


void gradient_3d(float *x, float *y, float *z, 
                 float *u, float *v, float *w, 
                 int i, int j, int k, 
                 int imax, int jmax, int kmax, 
                 float gradient_matrix[3][3])
{
  /// Find the gradient of a 3D vector.
  /// Returns a 3x3 gradient matrix.

  float u_xi, u_eta, u_zeta;
  float v_xi, v_eta, v_zeta;
  float w_xi, w_eta, w_zeta;
  float x_xi, x_eta, x_zeta;
  float y_xi, y_eta, y_zeta;
  float z_xi, z_eta, z_zeta;

  /// Using Finite difference scheme to calculate partial derivatives
  if (i == 0)
  {
    int index_i0 = get_index3(0,j,k,imax,jmax,kmax);
    int index_i1 = get_index3(1,j,k,imax,jmax,kmax);

    /// forward difference
    u_xi = *(u + index_i1) - *(u + index_i0);
    v_xi = *(v + index_i1) - *(v + index_i0);
    w_xi = *(w + index_i1) - *(w + index_i0);

    x_xi = *(x + index_i1) - *(x + index_i0);
    y_xi = *(y + index_i1) - *(y + index_i0);
    z_xi = *(z + index_i1) - *(z + index_i0);
  }
  else 
  { if (i == imax-1)
    {
      int index_imax_m1 = get_index3(imax-1,j,k,imax,jmax,kmax);
      int index_imax_m2 = get_index3(imax-2,j,k,imax,jmax,kmax);

      /// backward difference
      u_xi = *(u + index_imax_m1) - *(u + index_imax_m2);
      v_xi = *(v + index_imax_m1) - *(v + index_imax_m2);
      w_xi = *(w + index_imax_m1) - *(w + index_imax_m2);

      x_xi = *(x + index_imax_m1) - *(x + index_imax_m2);
      y_xi = *(y + index_imax_m1) - *(y + index_imax_m2);
      z_xi = *(z + index_imax_m1) - *(z + index_imax_m2);
    }
    else
    {
      int index_ip1 = get_index3(i+1,j,k,imax,jmax,kmax);
      int index_im1 = get_index3(i-1,j,k,imax,jmax,kmax);

      /// central difference
      u_xi = 0.5 * (*(u + index_ip1) - *(u + index_im1));
      v_xi = 0.5 * (*(v + index_ip1) - *(v + index_im1));
      w_xi = 0.5 * (*(w + index_ip1) - *(w + index_im1));

      x_xi = 0.5 * (*(x + index_ip1) - *(x + index_im1));
      y_xi = 0.5 * (*(y + index_ip1) - *(y + index_im1));
      z_xi = 0.5 * (*(z + index_ip1) - *(z + index_im1));
    }
  }
  

  if (j == 0) 
  {
    int index_j0 = get_index3(i,0,k,imax,jmax,kmax);
    int index_j1 = get_index3(i,1,k,imax,jmax,kmax);

    u_eta = *(u + index_j1) - *(u + index_j0);
    v_eta = *(v + index_j1) - *(v + index_j0);
    w_eta = *(w + index_j1) - *(w + index_j0);

    x_eta = *(x + index_j1) - *(x + index_j0);
    y_eta = *(y + index_j1) - *(y + index_j0);
    z_eta = *(z + index_j1) - *(z + index_j0);
  }
  else 
  {
    if (j == jmax-1)
    {
      int index_jmax_m1 = get_index3(i,jmax-1,k,imax,jmax,kmax);
      int index_jmax_m2 = get_index3(i,jmax-2,k,imax,jmax,kmax);

      u_eta = *(u + index_jmax_m1) - *(u + index_jmax_m2);
      v_eta = *(v + index_jmax_m1) - *(v + index_jmax_m2);
      w_eta = *(w + index_jmax_m1) - *(w + index_jmax_m2);

      x_eta = *(x + index_jmax_m1) - *(x + index_jmax_m2);
      y_eta = *(y + index_jmax_m1) - *(y + index_jmax_m2);
      z_eta = *(z + index_jmax_m1) - *(z + index_jmax_m2);
    }  
    else
    {
      int index_jp1 = get_index3(i,j+1,k,imax,jmax,kmax);
      int index_jm1 = get_index3(i,j-1,k,imax,jmax,kmax);

      u_eta = 0.5 * (*(u + index_jp1) - *(u + index_jm1));
      v_eta = 0.5 * (*(v + index_jp1) - *(v + index_jm1));
      w_eta = 0.5 * (*(w + index_jp1) - *(w + index_jm1));

      x_eta = 0.5 * (*(x + index_jp1) - *(x + index_jm1));
      y_eta = 0.5 * (*(y + index_jp1) - *(y + index_jm1));
      z_eta = 0.5 * (*(z + index_jp1) - *(z + index_jm1));
    }
  }

  if (k == 0)
  {
    int index_k0 = get_index3(i,j,0,imax,jmax,kmax);
    int index_k1 = get_index3(i,j,1,imax,jmax,kmax);

    u_zeta = *(u + index_k1) - *(u + index_k0);
    v_zeta = *(v + index_k1) - *(v + index_k0);
    w_zeta = *(w + index_k1) - *(w + index_k0);

    x_zeta = *(x + index_k1) - *(x + index_k0);
    y_zeta = *(y + index_k1) - *(y + index_k0);
    z_zeta = *(z + index_k1) - *(z + index_k0);
  }
  else
  {
    if (k == kmax-1) 
    {
      int index_kmax_m1 = get_index3(i,j,kmax-1,imax,jmax,kmax);
      int index_kmax_m2 = get_index3(i,j,kmax-2,imax,jmax,kmax);

      u_zeta = *(u + index_kmax_m1) - *(u + index_kmax_m2);
      v_zeta = *(v + index_kmax_m1) - *(v + index_kmax_m2);
      w_zeta = *(w + index_kmax_m1) - *(w + index_kmax_m2);

      x_zeta = *(x + index_kmax_m1) - *(x + index_kmax_m2);
      y_zeta = *(y + index_kmax_m1) - *(y + index_kmax_m2);
      z_zeta = *(z + index_kmax_m1) - *(z + index_kmax_m2);
    }
    else
    {
      int index_kp1 = get_index3(i,j,k+1,imax,jmax,kmax);
      int index_km1 = get_index3(i,j,k-1,imax,jmax,kmax);

      u_zeta = 0.5 * (*(u + index_kp1) - *(u + index_km1));
      v_zeta = 0.5 * (*(v + index_kp1) - *(v + index_km1));
      w_zeta = 0.5 * (*(w + index_kp1) - *(w + index_km1));

      x_zeta = 0.5 * (*(x + index_kp1) - *(x + index_km1));
      y_zeta = 0.5 * (*(y + index_kp1) - *(y + index_km1));
      z_zeta = 0.5 * (*(z + index_kp1) - *(z + index_km1));
    }
  }

  /// determinant of Jacobian
  float det =   x_xi * (y_eta*z_zeta-y_zeta*z_eta) 
              - x_eta * (y_xi*z_zeta-y_zeta*z_xi) 
              + x_zeta * (y_xi*z_eta-y_eta*z_xi);
  
  det = 1.0 / det;

  float xi_x = det*(y_eta*z_zeta - y_zeta*z_eta);
  float xi_y = det*(x_zeta*z_eta - x_eta*z_zeta);
  float xi_z = det*(x_eta*y_zeta - x_zeta*y_eta);

  float eta_x = det*(y_zeta*z_xi - y_xi*z_zeta);
  float eta_y = det*(x_xi*z_zeta - x_zeta*z_xi);
  float eta_z = det*(x_zeta*y_xi - x_xi*y_zeta);

  float zeta_x = det*(y_xi*z_eta - y_eta*z_xi);
  float zeta_y = det*(x_eta*z_xi - x_xi*z_eta);
  float zeta_z = det*(x_xi*y_eta - x_eta*y_xi);

  /// Assembling the velocity gradient tensor

  gradient_matrix[0][0] = u_xi*xi_x + u_eta*eta_x + u_zeta*zeta_x;   // dudx
  gradient_matrix[0][1] = u_xi*xi_y + u_eta*eta_y + u_zeta*zeta_y;   // dudy
  gradient_matrix[0][2] = u_xi*xi_z + u_eta*eta_z + u_zeta*zeta_z;   // dudz

  gradient_matrix[1][0] = v_xi*xi_x + v_eta*eta_x + v_zeta*zeta_x;    // dvdx
  gradient_matrix[1][1] = v_xi*xi_y + v_eta*eta_y + v_zeta*zeta_y;    // dvdy
  gradient_matrix[1][2] = v_xi*xi_z + v_eta*eta_z + v_zeta*zeta_z;    // dvdz

  gradient_matrix[2][0] = w_xi*xi_x + w_eta*eta_x + w_zeta*zeta_x;    // dwdx
  gradient_matrix[2][1] = w_xi*xi_y + w_eta*eta_y + w_zeta*zeta_y;    // dwdy
  gradient_matrix[2][2] = w_xi*xi_z + w_eta*eta_z + w_zeta*zeta_z;    // dwdz

}


void liutex(float velocity_gradient_tensor[3][3], float r[3])
{
  /* velocity_gradient_tensor = [  du/dx  du/dy  du/dz
                                   dv/dx  dv/dy  dv/dz
                                   dw/dx  dw/dy  dw/dz  ]
  */
  // r = liutex vector r.

  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;

  /// Extract the eigenvalues from the velocity gradient tensor. ///
  /// This part of the method finds the eigenvalues of the velocity gradient tensor.
  /// If you have better methods that do this for you, you can use them.

  float a[3][3];

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      a[i][j] = velocity_gradient_tensor[i][j];
    }
  }

  /// !----- Cardano's solution for the cubic equation -----
  /// !---------------------------------------------------------------------
  /// 
  /// ! cubic equation:
  /// ! x^3 + aa * x^2 + bb * x + cc = 0
  /// 
  /// ! coefficients of characteristic equation for the velocity gradient tensor.

  // Negative the Trace of velocity gradient.
  float aa = - a[0][0] - a[1][1] - a[2][2];

  float bb =   a[0][0] * a[1][1] - a[0][1] * a[1][0] 
             + a[0][0] * a[2][2] - a[0][2] * a[2][0] 
             + a[1][1] * a[2][2] - a[1][2] * a[2][1];

  // Negative the Determinate of velocity gradient.
  float cc = - a[0][0] * a[1][1] * a[2][2] 
             - a[0][1] * a[1][2] * a[2][0] 
             - a[0][2] * a[1][0] * a[2][1] 
             + a[0][2] * a[1][1] * a[2][0] 
             + a[0][1] * a[1][0] * a[2][2] 
             + a[0][0] * a[1][2] * a[2][1];

  /// delta is the discriminant of characteristic equation for the velocity graidient tensor.
  float qq = (aa * aa - 3.0 * bb) / 9.0;
  float rr = (2.0 * pow(aa, 3) - 9.0 * aa * bb + 27.0 * cc) / 54.0;

  float delta = pow(qq, 3) - pow(rr, 2);
    
  /// If the discriminant is less than 0 (delta < 0) then the velocity gradient tensor has one
  /// real eigenvalue and two complex conjugate eigenvalues and thus Liutex exists; 
  /// Else, the velocity gradient tensor has three real eigenvalues and Liutex is equal to 0.

  if (delta < 0.0)
  {
    float R = 0.0;

    float sign = 1.0;
    if (rr < 0)
    {
      sign = -1.0;
    }
    else if (rr > 0)
    {
      sign = 1.0;
    }
    else
    {
      sign = 0.0;
    }

    float aaaa = -sign * pow(fabs(rr) + sqrt(-delta), 1.0 / 3.0);

    float bbbb = 0.0;
    if (aaaa != 0.0)
    {
      bbbb = qq / aaaa;
    }

    /// Imaginary/complex conjugate eigenvalues of the velocity gradient tensor.
    /// imaginary_eigenvalue = lambda_c.
    /// lambda_cr + lambda_ci * i = real_part + imag_part * i.
    // float lambda_cr = -0.5 * (aaaa + bbbb) - aa / 3.0;
    float lambda_ci = 0.5 * sqrt(3.0) * (aaaa - bbbb);

    /// real_eig_val is the real eigenvalue of the velocity gradient tensor.
    float real_eig_val = aaaa + bbbb - aa / 3.0;

    //// Calculating the real right eigenvalue.
    float delta1 = (a[0][0] - real_eig_val) * (a[1][1] - real_eig_val) - a[1][0] * a[0][1];
    float delta2 = (a[1][1] - real_eig_val) * (a[2][2] - real_eig_val) - a[1][2] * a[2][1];
    float delta3 = (a[0][0] - real_eig_val) * (a[2][2] - real_eig_val) - a[0][2] * a[2][0];

    if (fabs(delta1) >= fabs(delta2) && fabs(delta1) >= fabs(delta3))
    {
      r[0] = (-(a[1][1] - real_eig_val) * a[0][2] + a[0][1] * a[1][2]) / delta1;
      r[1] = (a[1][0] * a[0][2] - (a[0][0] - real_eig_val) * a[1][2])  / delta1;
      r[2] = 1.0;
    }
    else if (fabs(delta2) >= fabs(delta1) && fabs(delta2) >= fabs(delta3))
    {
      r[0] = 1.0;
      r[1] = (-(a[2][2] - real_eig_val) * a[1][0] + a[1][2] * a[2][0]) / delta2;
      r[2] = (a[2][1] * a[1][0] - (a[1][1] - real_eig_val) * a[2][0])  / delta2;
    }
    else if (fabs(delta3) >= fabs(delta1) && fabs(delta3) >= fabs(delta2))
    {
      r[0] = (-(a[2][2] - real_eig_val) * a[0][1] + a[0][2] * a[2][1]) / delta3;
      r[1] = 1.0;
      r[2] = (a[2][0] * a[0][1] - (a[0][0] - real_eig_val) * a[2][1])  / delta3;
    }
    else
    {
      printf("ERROR: delta1, delta2, delta3\n");
      return;
    }

    //// Calculate the Liutex magnitude.

    /// Vorticity = w.
    float w[3];
    w[0] = a[2][1] - a[1][2];
    w[1] = a[0][2] - a[2][1];
    w[2] = a[1][0] - a[0][1];

    /// Dot product of vorticity and normalized real eigenvector of a.
    float w_dot_r = 0.0;
    for (int i = 0; i < 3; i++)
    {
      w_dot_r = w_dot_r + w[i] * r[i];
    }

    /// Direction Condition of Liutex magnitude.
    if (w_dot_r < 0.0)
    {
      r[0] = -1 * r[0];
      r[1] = -1 * r[1];
      r[2] = -1 * r[2];
    }

    /// Recalculate vorticity dot r.
    w_dot_r = 0.0;
    for (int i = 0; i < 3; i++)
    {
      w_dot_r = w_dot_r + w[i] * r[i];
    }

    /// Use explicit formula to calculate the Liutex magnitude R.
    R = w_dot_r - sqrt( pow(w_dot_r, 2) - 4.0 * pow(lambda_ci, 2) );

    r[0] = R * r[0];
    r[1] = R * r[1];
    r[2] = R * r[2];
  }

}

