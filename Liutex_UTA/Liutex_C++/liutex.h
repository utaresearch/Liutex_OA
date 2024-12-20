////////////////////////////////////////////////////////////////////////////////////
/// Liutex method C++ code.
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

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>


void liutex(std::vector<std::vector<double>> velocity_gradient_tensor, double& R, std::vector<double>& r)
{
    /* velocity_gradient_tensor = [  du/dx  du/dy  du/dz
                                    dv/dx  dv/dy  dv/dz
                                    dw/dx  dw/dy  dw/dz  ]
    */
    /// R = Liutex Magnitude
    /// r = Liutex Vector

    /// Extract the eigenvalues from the velocity gradient tensor. ///
    /// This part of the method finds the eigenvalues of the velocity gradient tensor.
    /// If you have methods that do this for you, please use them.

    const std::vector<std::vector<double>> a = velocity_gradient_tensor;

    /// ! Cubic Formula
    /// ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
    /// ! 5.6 Quadratic and Cubic Equations
    /// ! Page 179
    /// !---------------------------------------------------------------------
    /// 
    /// ! cubic equation
    /// ! x**3 + aa * x**2 + bb * x + cc = 0
    /// 
    /// ! coefficients of characteristic equation for the velocity gradient tensor.

    double aa = -( a[0][0] + a[1][1] + a[2][2] );

    /// If you can find a method for multiplying matrices, you can use it.
    /// Squaring the velocity gradient tensor (i.e., tt = a^2 = a*a).
    std::vector<std::vector<double>> tt(3, std::vector<double>(3, 0.0));

    for (uint8_t i = 0; i < 3; i++)
    {
        for (uint8_t j = 0; j < 3; j++)
        {
            for (uint8_t k = 0; k < 3; k++)
            {
                tt[i][j] += a[i][k] * a[k][j];
            }
        }
    }

    double bb = -0.5 * ( tt[0][0] + tt[1][1] + tt[2][2] - pow(a[0][0] + a[1][1] + a[2][2], 2) );

    double cc = -( a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
                - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
                + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]) );

    /// delta is the discriminant of characteristic equation for the velocity graidient tensor.
    double delta = 18.0*aa*bb*cc - 4.0*pow(aa,3) * cc + aa*aa * bb*bb - 4.0*pow(bb,3) - 27.0*cc*cc;
    delta = delta / 108.0;
    
    /// If the discriminant is less than 0 (delta < 0) then the velocity gradient tensor has one
    /// real eigenvalue and two complex conjugate eigenvalues and thus Liutex exists. 
    /// Else, the velocity gradient tensor has three real eigenvalues and Liutex is equal to 0.
    if (delta < 0.0)
    {
        delta = -delta;

        double qq = (aa*aa - 3.0*bb) / 9.0;
        double rr = (2.0*pow(aa,3) - 9.0*aa*bb + 27.0*cc) / 54.0;

        double sign;
        if (std::signbit(rr))
        {
            sign = -1.0;
        }
        else
        {
            sign = 1.0;
        }
        double aaaa = -sign * pow(abs(rr) + sqrt(delta), 1.0/3.0);

        double bbbb;
        if (aaaa == 0.0)
        {
            bbbb = 0.0;
        }
        else
        {
            bbbb = qq / aaaa;
        }

        /// eig1c, eig2c = the complex conjugate eigenvalues of the velocity gradient tensor.
        /// eig3r = the real eigenvalue of the velocity graident tensor.
        std::complex<double> eig1c( -0.5*(aaaa+bbbb) - aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb) );
        std::complex<double> eig2c( eig1c.real(), -eig1c.imag() );
        double eig3r = aaaa + bbbb - aa/3.0;

        /// Calculating the real right eigenvalue.
        double delta1 = (a[0][0] - eig3r) * (a[1][1] - eig3r) - a[1][0]*a[0][1];
        double delta2 = (a[1][1] - eig3r) * (a[2][2] - eig3r) - a[1][2]*a[2][1];
        double delta3 = (a[0][0] - eig3r) * (a[2][2] - eig3r) - a[0][2]*a[2][0];

        if (abs(delta1) >= abs(delta2) && abs(delta1) >= abs(delta3))
        {
            r[0] = (-(a[1][1]-eig3r)*a[0][2] +         a[0][1]*a[1][2]) / delta1;
            r[1] = (         a[1][0]*a[0][2] - (a[0][0]-eig3r)*a[1][2]) / delta1;
            r[2] = 1.0;
        }
        else if (abs(delta2) >= abs(delta1) && abs(delta2) >= abs(delta3))
        {
            r[0] = 1.0;
            r[1] = (-(a[2][2]-eig3r)*a[1][0] +         a[1][2]*a[2][0])/delta2;
            r[2] = (         a[2][1]*a[1][0] - (a[1][1]-eig3r)*a[2][0])/delta2;
        }
        else if (abs(delta3) >= abs(delta1) && abs(delta3) >= abs(delta2))
        {
            r[0] = (-(a[2][2]-eig3r)*a[0][1] +         a[0][2]*a[2][1])/delta3;
            r[1] = 1.0;
            r[2] = (         a[2][0]*a[0][1] - (a[0][0]-eig3r)*a[2][1])/delta3;
        }
        else
        {
            std::cout << "ERROR: delta1, delta2, delta3\n";
            return;
        }

        double r_norm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

        r[0] = r[0] / r_norm;
        r[1] = r[1] / r_norm;
        r[2] = r[2] / r_norm;

        /// Calculate rotation matrix r which rotates unit vector u to unit vector v.

        /// u = z0 = 0, 0, 1; 
        /// v = r = liutex_vec
        /// r = qqq= transformation matrix

        std::vector<double> z0 = { 0.0, 0.0, 1.0 };
        std::vector<double> temp_vec(3);
        std::vector<std::vector<double>> transformation_matrix(3, std::vector<double> (3, 0.0));

        double eps = 1.0e-10;

        temp_vec[0] = z0[1]*r[2] - z0[2]*r[1];
        temp_vec[1] = z0[2]*r[0] - z0[0]*r[2];
        temp_vec[2] = z0[0]*r[1] - z0[1]*r[0];

        aa = sqrt( temp_vec[0]*temp_vec[0] + temp_vec[1]*temp_vec[1] + temp_vec[2]*temp_vec[2] );
        
        if (aa < eps)
        {
            transformation_matrix[0][0] = 1.0;
            transformation_matrix[1][1] = 1.0;
            transformation_matrix[2][2] = 1.0;
        }
        else
        {
            temp_vec[0] = temp_vec[0] / aa;
            temp_vec[1] = temp_vec[1] / aa;
            temp_vec[2] = temp_vec[2] / aa;

            double t = z0[0]*r[0] + z0[1]*r[1] + z0[2]*r[2];
            
            if (t > 1.0)
            {
                t = 1.0;
            }

            if (t < -1.0)
            {
                t = -1.0;
            }
            
            double alpha = acos(t);
        
            double c = cos(alpha);
            double s = sin(alpha);
        
            transformation_matrix[0][0] = temp_vec[0] * temp_vec[0] * (1.0 - c) + c;
            transformation_matrix[0][1] = temp_vec[0] * temp_vec[1] * (1.0 - c) - temp_vec[2] * s;
            transformation_matrix[0][2] = temp_vec[0] * temp_vec[2] * (1.0 - c) + temp_vec[1] * s;
        
            transformation_matrix[1][0] = temp_vec[1] * temp_vec[0] * (1.0 - c) + temp_vec[2] * s;
            transformation_matrix[1][1] = temp_vec[1] * temp_vec[1] * (1.0 - c) + c;
            transformation_matrix[1][2] = temp_vec[1] * temp_vec[2] * (1.0 - c) - temp_vec[0] * s;
        
            transformation_matrix[2][0] = temp_vec[2] * temp_vec[0] * (1.0 - c) - temp_vec[1] * s;
            transformation_matrix[2][1] = temp_vec[2] * temp_vec[1] * (1.0 - c) + temp_vec[0] * s;
            transformation_matrix[2][2] = temp_vec[2] * temp_vec[2] * (1.0 - c) + c;
        }

        /// If you can find a method for multiplying matrices, you can use it.
        /// vg = transpose(transformation_matrix) * velocity_gradient_tensor * transformation_matrix.
        
        /// Transpose the transformation matrix.
        std::vector<std::vector<double>> transpose_transformation_matrix(3, std::vector<double> (3, 0.0));
        
        for (uint8_t i = 0; i < 3 ; i++)
        {
            for (uint8_t j = 0; j < 3 ; j++)
            {
                transpose_transformation_matrix[i][j] = transformation_matrix[j][i];
            }

        }

        /// matrix_product_1 = transpose_transformation_matrix * velocity_gradient_tensor.
        std::vector<std::vector<double>> matrix_product_1(3, std::vector<double>(3, 0.0));

        for (uint8_t i = 0; i < 3; i++)
        {
            for (uint8_t j = 0; j < 3; j++)
            {
                matrix_product_1[i][j] = 0.0;

                for (uint8_t k = 0; k < 3; k++)
                {
                    matrix_product_1[i][j] += transpose_transformation_matrix[i][k] * a[k][j];
                }
            }
        }

        /// vg = velocity_gradient_tensor * transformation_matrix.
        std::vector<std::vector<double>> vg(3, std::vector<double>(3, 0.0));

        for (uint8_t i = 0; i < 3; i++)
        {
            for (uint8_t j = 0; j < 3; j++)
            {
                vg[i][j] = 0.0;

                for (uint8_t k = 0; k < 3; k++)
                {
                    vg[i][j] += a[i][k] * transformation_matrix[k][j];
                }
            }
        }

        double alpha = sqrt( pow(vg[1][1] - vg[0][0], 2) + pow(vg[1][0] + vg[0][1], 2) );
        double beta  = (vg[1][0] - vg[0][1]);

        if (beta*beta > alpha*alpha)
        {
            if(beta > 0.0)
            {
                R = 2.0 * (beta - alpha);
                r[0] = R * r[0];
                r[1] = R * r[1];
                r[2] = R * r[2];
            }
            else
            {
                R = 2.0 * (beta + alpha);
                r[0] = R * r[0];
                r[1] = R * r[1];
                r[2] = R * r[2];
            }
        }
        else
        {
            r[0] = 0.0;
            r[1] = 0.0;
            r[2] = 0.0;
        }

        R = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
    }
    else
    {
        /// All eigenvalue are real so Liutex = 0.
        R = 0.0;
        r[0] = 0.0;
        r[1] = 0.0;
        r[2] = 0.0;
    }



    return;

}

