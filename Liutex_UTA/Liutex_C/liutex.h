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


void liutex(double velocity_gradient_tensor[3][3], double r[3])
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

  double a[3][3] = { 0.0 };

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
  double aa = - a[0][0] - a[1][1] - a[2][2];

  double bb =   a[0][0] * a[1][1] - a[0][1] * a[1][0] 
             + a[0][0] * a[2][2] - a[0][2] * a[2][0] 
             + a[1][1] * a[2][2] - a[1][2] * a[2][1];

  // Negative the Determinate of velocity gradient.
  double cc = - a[0][0] * a[1][1] * a[2][2] 
             - a[0][1] * a[1][2] * a[2][0] 
             - a[0][2] * a[1][0] * a[2][1] 
             + a[0][2] * a[1][1] * a[2][0] 
             + a[0][1] * a[1][0] * a[2][2] 
             + a[0][0] * a[1][2] * a[2][1];

  /// delta is the discriminant of characteristic equation for the velocity graidient tensor.
  double qq = (aa * aa - 3.0 * bb) / 9.0;
  double rr = (2.0 * pow(aa, 3) - 9.0 * aa * bb + 27.0 * cc) / 54.0;

  double delta = pow(qq, 3) - pow(rr, 2);
    
  /// If the discriminant is less than 0 (delta < 0) then the velocity gradient tensor has one
  /// real eigenvalue and two complex conjugate eigenvalues and thus Liutex exists; 
  /// Else, the velocity gradient tensor has three real eigenvalues and Liutex is equal to 0.

  if (delta < 0.0)
  {
    double R = 0.0;

    double sign = 1.0;
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

    double aaaa = -sign * pow(abs(rr) + sqrt(-delta), 1.0 / 3.0);

    double bbbb = 0.0;
    if (aaaa != 0.0)
    {
      bbbb = qq / aaaa;
    }

    /// Imaginary/complex conjugate eigenvalues of the velocity gradient tensor.
    /// imaginary_eigenvalue = lambda_c.
    /// lambda_cr + lambda_ci * i = real_part + imag_part * i.
    // double lambda_cr = -0.5 * (aaaa + bbbb) - aa / 3.0;
    double lambda_ci = 0.5 * sqrt(3.0) * (aaaa - bbbb);

    /// real_eig_val is the real eigenvalue of the velocity gradient tensor.
    double real_eig_val = aaaa + bbbb - aa / 3.0;

    //// Calculating the real right eigenvalue.
    double delta1 = (a[0][0] - real_eig_val) * (a[1][1] - real_eig_val) - a[1][0] * a[0][1];
    double delta2 = (a[1][1] - real_eig_val) * (a[2][2] - real_eig_val) - a[1][2] * a[2][1];
    double delta3 = (a[0][0] - real_eig_val) * (a[2][2] - real_eig_val) - a[0][2] * a[2][0];

    if (abs(delta1) >= abs(delta2) && abs(delta1) >= abs(delta3))
    {
      r[0] = (-(a[1][1] - real_eig_val) * a[0][2] + a[0][1] * a[1][2]) / delta1;
      r[1] = (a[1][0] * a[0][2] - (a[0][0] - real_eig_val) * a[1][2])  / delta1;
      r[2] = 1.0;
    }
    else if (abs(delta2) >= abs(delta1) && abs(delta2) >= abs(delta3))
    {
      r[0] = 1.0;
      r[1] = (-(a[2][2] - real_eig_val) * a[1][0] + a[1][2] * a[2][0]) / delta2;
      r[2] = (a[2][1] * a[1][0] - (a[1][1] - real_eig_val) * a[2][0])  / delta2;
    }
    else if (abs(delta3) >= abs(delta1) && abs(delta3) >= abs(delta2))
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
    double w[3] = { 0.0 };
    w[0] = a[2][1] - a[1][2];
    w[1] = a[0][2] - a[2][1];
    w[2] = a[1][0] - a[0][1];

    /// Dot product of vorticity and normalized real eigenvector of a.
    double w_dot_r = 0.0;
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

