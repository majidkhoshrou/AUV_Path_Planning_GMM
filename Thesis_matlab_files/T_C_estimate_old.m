function [C, T, D] = T_C_estimate(in)
%#codegen
% x can anything
% y [y_left_d y_right_d]
z = in(1);
x = in(2);
y = in(3);
y_left_d = in(4);
y_right_d = in(5);
D = z * 11.55 ;
y_original = y ;
y = (y - y_left_d)/(y_right_d - y_left_d) ;
% pi = 3.1415 ;
...........................................................................
% popup4 regression model
% temperature
% Coefficients (with 95% confidence bounds):
       a =       19.91  ; % (19.9, 19.91)
       b =    -0.02081  ; % (-0.02687, -0.01475)
       c =      0.2248  ; % (0.2138, 0.2358)
       m =       2.879  ; % (2.695, 3.063)
       w =       2.615  ; % (2.468, 2.763)
T_right = a + b*sin(m*pi*y*z) + c*exp(-(w*z)^2) ;
...
    % Coefficients (with 95% confidence bounds):
       a =       4.938  ; % (-2.58e+04, 2.581e+04)
       b =     -0.0162  ; % (-0.0186, -0.01379)
       c =     -0.0148  ; % (-2.58e+04, 2.58e+04)
       m =      0.8842  ; % (0.8181, 0.9504)
       w =   -0.006801  ; % (-5928, 5928)
C_right = a + b*sin(m*pi*y*z) + c*exp(-(w*z)^2) ;

...........................................................................
% popup2 regression model
        % Coefficients (with 95% confidence bounds):
% Temperature
       a =       19.2  ; % (19.24, 19.46)   19.35
       b =     0.01291  ; % (0.008128, 0.0177)
       c =      0.4754  ; % (0.3695, 0.5813)
       m =       2.851  ; % (2.641, 3.061)
       w =     -0.8884  ; % (-1.041, -0.7359)
T_left = a + b*sin(m*pi*y*z) + c*exp(-(w*z)^2) ;

...
  % Coefficients (with 95% confidence bounds):
       a =       4.7  ; % (4.763, 4.864)   4.813
       b =   -0.003453  ; % (-0.004459, -0.002447)
       c =     0.08822  ; % (0.03802, 0.1384)
       m =       1.305  ; % (1.176, 1.434)
       w =      0.6251  ; % (0.4072, 0.8429)
C_left = a + b*sin(m*pi*y*z) + c*exp(-(w*z)^2) ;



                 if y_original <= y_left_d
                     C = C_left ;
                     T = T_left ;
                 elseif y_original >= y_right_d
                     C = C_right ;
                     T = T_right ;
                 else                     
                     C = ((y_right_d - y_original)* C_left + (y_original - y_left_d)* C_right)/(y_right_d-y_left_d) ;
                     T = ((y_right_d - y_original)* T_left + (y_original - y_left_d)* T_right)/(y_right_d-y_left_d) ;
                 end
                 

                 
                 