clear all
close all
clc

% System Parameters
n = 4;  % Number of Vessels
T_f = 1;  % Formation Control Frequency
T_bc =  1;   % Broadcast Intervall
T_est = .1;  % Esimtator sampling time

%Initial position [P0 psi0]
P01 = [20   20];
% P02 = [0   -10];
P02 = [40   00];
P03 = [20  -20];
P04 = [150 40];
psi01 = 0*pi/180;
psi02 = 0*pi/180;
psi03 = 0*pi/180;
psi04 = 0*pi/180;
loc1 = [0 0]';
loc2 = [0 -5]';
loc3 = [0 5]';
loc4 = [0 0]';
D = [2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 0];
A = [0 1 1 0; 1 0 1 0; 1 1 0 0; 0 0 0 0];
ksync = .01;

%Initial speed [u0 v0 r0]
U01 = [.5 0 0]';
U02 = [.5 0 0]';
U03 = [.5 0 0]';

%Sensor noise power
noiseon = 1;
npsi = 0.00005;
nu = .05;
nv = 0.005;
nr = 0.005;
% nsonar = 0.05;
nsonar = .1;
nWvxy = .5;

%Estimator
sN = 0.1*[ 1.0753    3.6678   -4.5177    1.7243    0.6375   -2.6154];
% SN(1:2) = [2 2];
State_est01 = [P01+[sN(1),sN(2)] psi01];
State_est02 = [P02+[sN(3),sN(4)] psi02];
State_est03 = [P03+[sN(5),sN(6)] psi03];
Sigma01 = 15*eye(3);
Sigma02 = Sigma01;
Sigma03 = Sigma01;

% Plant noise covariance matrices
Q1 = diag([nu nv npsi]);
Q2 = diag([nu nv npsi]);
Q3 = diag([nu nv npsi]);

% Transmission laws
% TrsmtLaw = 4;
%   TrsmtLaw = [1  2  3];
%   TrsmtLaw = [4 1  3];
% TrsmtLaw = [4 1  2  3];
TrsmtLaw = [4  1  4  2  4  3];
% T_bc =  100000;

% Pathdata for all vehicles
%Vehicle 1          gamWP   [yaw/0]   [ P0/R]    v_L  Radius
Pathdata(:,1,1) = [  0      0           P01      .5    0]';
Pathdata(:,2,1) = [500     -pi/2     270  25     .5    5]';
Pathdata(:,3,1) = [531.4    pi       270  30     .5    0]';
Pathdata(:,4,1) = [611      pi       230  30     .5    0]';
Pathdata(:,5,1) = [2998     pi       230  30     .5    0]';
Pathdata(:,6,1) = [2999     pi       230  30     .5    0]';
% %Vehicle 2          gamWP   [theta]                 [ P0]     v_L
% Pathdata(:,1,2) = [  0      0           P02      .5    0]';
% Pathdata(:,2,2) = [ 20     -pi/2      10   0     .5   10]';
% Pathdata(:,3,2) = [2996     pi        00   0     .5    0]';
% Pathdata(:,4,2) = [2997     pi        00   0     0     0]';
% Pathdata(:,5,2) = [2998     pi        00   0     0     0]';
% Pathdata(:,6,2) = [2999     pi        00   0     0     0]';
%Vehicle 2          gamWP   [theta]                 [ P0]    v_L
Pathdata(:,1,2) = [  0      0           P02      .5    0]';
Pathdata(:,2,2) = [500     -pi/2     290   5     .5    5]';
Pathdata(:,3,2) = [531.4    pi       290  10     .5    0]';
Pathdata(:,4,2) = [611      pi       250  10     .5    0]';
Pathdata(:,5,2) = [2998     pi       250  10     .5    0]';
Pathdata(:,6,2) = [2999     pi       250  10     .5    0]';
%Vehicle 3          gamWP   [theta]                 [ P0]    v_L
Pathdata(:,1,3) = [  0      0           P03      .5    0]';
Pathdata(:,2,3) = [500     -pi/2     270 -15     .5    5]';
Pathdata(:,3,3) = [531.4    pi       270 -10     .5    0]';
Pathdata(:,4,3) = [611      pi       230 -10     .5    0]';
Pathdata(:,5,3) = [2998     pi       230 -10     .5    0]';
Pathdata(:,6,3) = [2999     pi       230 -10     .5    0]';
% Beacon 1         gamWP   [theta]              [ P0]    v_L
Pathdata(:,1,4) = [   0     0           P04       0    0]';
Pathdata(:,2,4) = [  20     0           P04       0    0]';
Pathdata(:,3,4) = [  30     0           P04       0    0]';
Pathdata(:,4,4) = [  40     0           P04       0    0]';
Pathdata(:,5,4) = [  50     0           P04       0    0]';
Pathdata(:,6,4) = [ 2999    0           P04       0    0];

%Syrene 
rho = 1026;     %kg/m^3
m = 4000;       %kg
V = 3.998;      %m^3
I_z = 2657.5;   %kg*m^2

%Added mass
Xm_u = -0.07*rho*V;
Xm_v = -0.077*rho*V;
N = -95;% -0.085*rho*V^(5/3);
m_u = m - Xm_u;
m_v = m - Xm_v;
m_uv = Xm_v - Xm_u;
J = I_z - N;

%Viscous forces
Xu = -360;      %kg/s
Xuu = -805;     %kg/m
Yv = -420;      %kg/s
Yvv = -1930;    %kg/m
Nr = -110;      %kg*m/s
Nrr = -555;     %kg*m

%Propulsion
Fmax = 1000;    %N (1500, -1000)
Fmin = -750;    %N
l = 1;          %m

%Ocean Current (with normal can be up to 1 m/s)
Vc = .3;
% Vc = 0;
psi_c = 45*pi/180;
nVc = 0*.005;
npsi_c = 0*.00005;

%Inner Loop Control
k_u = .5;
k_r = .5;

%Outer Loop Control (Normal: .3 .15 1.5 testare kx con waypoint - BS: 0 0 1.5 kw = [10 0; 0 1];)
delta = .5;
kx = .3;
ky = .15;
kz = 1.5;
near = .7;
kw = [10 0; 0 1];
Ksync = diag(ksync*ones(1,n));
Kz = diag(kz*ones(1,n));

%Prefilter
tau_1 = .1;
tau_2 = 11;
T1 = 10;
T2 = 2;

% Coordination Controller
k_xi = 1;