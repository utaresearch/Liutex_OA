%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liutex method MATLAB code.
% This code calculates liutex when given the veloctiy gradient tensor 
% for 3D flow fields.
% 
% Author:  Oscar Alvarez
% 
% email: oscar.alvarez@uta.edu
% University of Texas at Arlington
% Department of Mathematics
% Center for Numerical Simulation and Modeling (CNSM)
% Arlington, Texas, United States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, r] = liutex(velocity_gradient_tensor)
% This function returns the Liutex Magnitude R and
% the Liutex vector r when given the velocity gradient
% tensor (3x3 matrix) for a 3-dimensional flow field.

    A = velocity_gradient_tensor;
    
    % Extract the eigenvalues from the velocity gradient tensor.
    [eig_vec, eig_val] = eig(A);
    
    % Check if the eigenvalues of the velocity gradient tensor
    % are all real. If they are real, then Liutex is = 0.
    if isreal(eig_val)
        R = 0;
        r = [0, 0, 0];
        return;
    end
    
    % If there exist complex eigenvalues, then continue to
    % calculate Liutex.

    % Determine which eigenvalues/eigenvectors are real and 
    % which are complex.
    for i = 1:3
        if (imag(eig_val(i,i)) == 0)
            real_index = i;
            imaginary_index = mod(i, 3) + 1;
        end
    end

    % Save the real eigenvector of the velocity gradient tensor.
    real_eig_vec = real( eig_vec(:,real_index) );

    % Save the imaginary part of the complex eigenvalues of the 
    % velocity gradient tensor.
    lambda_ci = imag( eig_val(imaginary_index,imaginary_index) );
    
    % Calculate vorticity using the velocity gradient tensor.
    vorticity(1) = A(3,2) - A(2,3);
    vorticity(2) = A(1,3) - A(3,1);
    vorticity(3) = A(2,1) - A(1,2);
    
    % Check first condition of Liutex and set r.
    if (dot(vorticity, real_eig_vec) < 0)
        r = -real_eig_vec;
    else
        r = real_eig_vec;
    end
    
    % Now use explicit formula to calculate R.
    w_dot_r = dot(vorticity, r);
    R = w_dot_r - sqrt( w_dot_r^2 - 4*lambda_ci^2 );
    
    % Apply R to r and finalize the Liutex vector r.
    r = R * r;
    
    return;
end