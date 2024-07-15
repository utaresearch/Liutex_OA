////////////////////////////////////////////////////////////////////////////////////
/// Liutex example C code.
/// This code shows how to use the liutex method located in the "liutex.h" file.
/// 
/// Author:  Oscar Alvarez
/// 
/// email: oscar.alvarez@uta.edu
/// University of Texas at Arlington
/// Department of Mathematics
/// Center for Numerical Simulation and Modeling (CNSM)
/// Arlington, Texas, United States
////////////////////////////////////////////////////////////////////////////////////

#include "liutex.h"

#include <stdio.h>

int main()
{
    printf("Math Department - University of Texas Arlington \n");
    printf("\nTesting Liutex Function with sample velocity gradient tensors.\n");
    
    /// real_eigvals_vgt := velocity gradient tensor with all real eigenvalues.
    /// imag_eigvals_vgt := velocity gradient tensor with imaginary eigenvalues.
    double real_eigvals_vg[3][3] = { {-2.0, -4.0, 2.0}, 
                                     {-2.0,  1.0, 2.0},
                                     { 4.0,  2.0, 5.0} };

    double imag_eigvals_vg[3][3] = { {4.0, -3.0,  7.0}, 
                                     {3.0,  4.0,  0.0},
                                     {5.0, 10.0, 10.0} };

    /// Initializing Liutex vector (r).
    double r_real[3];
    double r_imag[3];
    
    liutex(real_eigvals_vg, r_real);

    /// Calculate Liutex magnitude.
    double R_real = sqrt(r_real[0] * r_real[0] + r_real[1] * r_real[1] + r_real[2] * r_real[2]);

    /// Printing Results:
    printf("\n Velocity gradient tensor with all real eigenvalues \n");
    printf("Liutex magnitude R: %f \n", R_real);
    printf("Liutex vector r:  \n");
    for (int i = 0; i < 3; i++)
    {
        printf("%f \n", r_real[i]);
    }

    liutex(imag_eigvals_vg, r_imag);

    /// Calculate Liutex magnitude.
    double R_imag = sqrt(r_imag[0] * r_imag[0] + r_imag[1] * r_imag[1] + r_imag[2] * r_imag[2]);

    /// Printing Results:
    printf("\n\n Velocity gradient tensor with complex conjugate eigenvalues \n");
    printf("Liutex magnitude R: %f \n", R_imag);
    printf("Liutex vector r:  \n");
    for (int i = 0; i < 3; i++)
    {
        printf("%f \n", r_imag[i]);
    }

    printf("\n");

    return 0;

}