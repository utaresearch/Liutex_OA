////////////////////////////////////////////////////////////////////////////////////
/// Liutex example C++ code.
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

#include <iostream>
#include <vector>

int main()
{
    std::cout << "Math Department - University of Texas Arlington \n";
    std::cout << "\nTesting Liutex Function \n";
    
    /// real_eigvals_vgt := velocity gradient tensor with all real eigenvalues.
    /// imag_eigvals_vgt := velocity gradient tensor with imaginary eigenvalues.
    std::vector<std::vector<double>> real_eigvals_vgt;
    std::vector<std::vector<double>> imag_eigvals_vgt;

    real_eigvals_vgt = { {-2.0, -4.0, 2.0}, {-2.0, 1.0, 2.0}, {4.0, 2.0, 5.0} };
    imag_eigvals_vgt = { {4.0, -3.0, 7.0}, {3.0, 4.0, 0.0}, {5.0, 10.0, 10.0} };

    /// Initializing Liutex magnitude (R) and Liutex vector (r).
    double R_real, R_imag;
    std::vector<double> r_real(3), r_imag(3);
    
    liutex(real_eigvals_vgt, R_real, r_real);
    liutex(imag_eigvals_vgt, R_imag, r_imag);

    /// Printing Results:
    std::cout << "\n[ Velocity gradient tensor with all real eigenvalues ]\n";
    std::cout << "Liutex magnitude R: " << R_real << "\n";
    std::cout << "Liutex vector r:  ";
    for (auto x : r_real)
    {
        std::cout << x << "  ";
    }

    std::cout << "\n\n[ Velocity gradient tensor with complex conjugate eigenvalues ]\n";
    std::cout << "Liutex magnitude R: " << R_imag << "\n";
    std::cout << "Liutex vector r:  ";
    for (auto x : r_imag)
    {
        std::cout << x << "  ";
    }

    std::cout << "\n";

    return 0;

}