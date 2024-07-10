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
    printf("\nTesting Liutex Function \n");
    
    /// real_eigvals_vgt := velocity gradient tensor with all real eigenvalues.
    /// imag_eigvals_vgt := velocity gradient tensor with imaginary eigenvalues.
    double real_eigvals_vgt[3][3] = { {-2.0, -4.0, 2.0}, {-2.0, 1.0, 2.0}, {4.0, 2.0, 5.0} };
    double imag_eigvals_vgt[3][3] = { {4.0, -3.0, 7.0}, {3.0, 4.0, 0.0}, {5.0, 10.0, 10.0} };

    /// Initializing Liutex magnitude (R) and Liutex vector (r).
    double R_real, R_imag;
    double r_real[3], r_imag[3];
    
    liutex(real_eigvals_vgt, &R_real, &r_real);
    liutex(imag_eigvals_vgt, &R_imag, &r_imag);

    /// Printing Results:
    printf("\n\[ Velocity gradient tensor with all real eigenvalues \]\n");
    printf("Liutex magnitude R: ", R_real, "\n");
    printf("Liutex vector r:  ");
    for (int i = 0; i < sizeof(r_real); i++)
    {
        printf(" ", r_real[i], " ");
    }

    printf("\n\n[ Velocity gradient tensor with complex conjugate eigenvalues ]\n");
    printf("Liutex magnitude R: ", R_imag, "\n");
    printf("Liutex vector r:  ");
    for (int i = 0; i < sizeof(r_imag); i++)
    {
        printf(" ", r_imag[i], "  ");
    }

    printf("\n");

    return 0;

}