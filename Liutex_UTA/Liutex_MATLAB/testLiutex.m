%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liutex example MATLAB code.
% This code shows how to use the liutex function located in the "liutex.m" file.
% 
% Author:  Oscar Alvarez
% 
% email: oscar.alvarez@uta.edu
% University of Texas at Arlington
% Department of Mathematics
% Center for Numerical Simulation and Modeling (CNSM)
% Arlington, Texas, United States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; clc;

% Velocity gradient tensor with only real eigenvalues.
real_velocity_gradient_tensor = [-2 -4 2 ; -2 1 2 ; 4 2 5];

% Velocity gradient tensor with complex eigenvalues.
imag_velocity_gradient_tensor = [4.0 -3.0 7.0 ; 3.0 4.0 0.0 ; 5.0 10.0 10.0];


fprintf("Results:\n")

% Test velocity gradient tensor with real eigenvalues.
[R, r] = liutex(real_velocity_gradient_tensor);
fprintf("Real eigenvalues: \n");
fprintf("R= %g \n", R);
fprintf("r= %g \n", r);

fprintf("\n");

% Test velocity gradient tensor with complex eigenvalues.
[R, r] = liutex(imag_velocity_gradient_tensor);
fprintf("Imaginary eigenvalues: \n")
fprintf("R= %g \n", R);
fprintf("r= %g \n",r);
