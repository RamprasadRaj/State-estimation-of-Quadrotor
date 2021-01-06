% X config
clear
clc
%
mass = 0.25; % mass of quadcopter
iX = 0.03; % kgm^2
iY = 0.03; % kgm^2
iZ = 0.03; % kgm^2
% thrustCoeff = 0.1;
% rho = 1.225;
% rotorRadius = 0.035;
% rotorArea = pi*rotorRadius^2;
% thrustFactor = thrustCoeff*rho*rotorArea*rotorRadius^2;
thrustFactor = 0.1;
dragFactor = 0.1;
armLength = 0.1;
%% Initial conditions
%
x0 = 0.0;
y0 = 0.0;
z0 = 0.0;
xDot0 = 0.0;
yDot0 = 0.0;
zDot0 = 0.0;
phi0 = 0.0;
theta0 = 0.0;
psi0 = 0.0;
phiDot0 = 0.0;
thetaDot0 = 0.0;
psiDot0 = 0.0;
%
%% Input
% The inputs are propeller speeds and are given such that the quadcopter
% pitches forward and move along x-axis while climbing up.
om1 = 2.8;
om2 = 3.2;
om3 = 3.2;
om4 = 2.8;
%
%% Extended Kalman filter
%
EKFquadcopter
%% Unscented Kalman Filter
UKFquadcopter