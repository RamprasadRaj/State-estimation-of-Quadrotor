%% Extended Kalman Filter
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
Ts = 0.025;
totT = 2.5;
sigmaGPS = 0.5;
sigmaGY = 0.1;
e = 1e-6;
x_k0 = [x0;y0;z0;phi0;theta0;psi0;...
    xDot0;yDot0;zDot0;phiDot0;thetaDot0;psiDot0];
count = 1;
t = Ts:Ts:totT;
P_k1 = eye(12);
Q = diag(e*ones(1,12));
R = eye(6);
R(1,1) = sigmaGPS*R(1,1);
R(2,2) = sigmaGPS*R(2,2);
R(3,3) = sigmaGPS*R(3,3);
R(4,4) = sigmaGY*R(4,4);
R(5,5) = sigmaGY*R(5,5);
R(6,6) = sigmaGY*R(6,6);
H = [1 0 0 0 0 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0 0 0 0 0;...
    0 0 1 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 1 0 0;...
    0 0 0 0 0 0 0 0 0 0 1 0;...
    0 0 0 0 0 0 0 0 0 0 0 1];
simout = sim('EENGquadcopter_control','StartTime',num2str(0),'StopTime',...
        num2str(Ts),'OutputOption','SpecifiedOutputTimes','OutputTimes',...
        num2str(Ts));
input = simout.Input.Data(end,:)';
del_x_kp1(:,count) = simout.simout.Data(end,:)';
yMeasurement(:,count) = simout.sensorMeasure.Data(end,:)';
yMeasurement(:,count) = yMeasurement(:,count) + [sigmaGPS*randn(3,1);...
    sigmaGY*randn(3,1)];
%
deltaF = zeros(12); % Derivative of the non-linear function matrix
%
F = zeros(12,1); % Non-linear function matrix
%
deltaF(1,7) = 1.0; % xDot
deltaF(2,8) = 1.0; % yDot
deltaF(3,9) = 1.0; % zDot
deltaF(4,10) = 1.0; % phiDot
deltaF(5,11) = 1.0; % thetaDot
deltaF(6,12) = 1.0; % psiDot
%
deltaF(7,4) = (input(1)/mass)*((-sin(phi0)*sin(theta0)*cos(psi0))+(cos(phi0)*sin(psi0)));
deltaF(7,5) = (input(1)/mass)*(cos(phi0)*cos(theta0)*cos(psi0));
deltaF(7,6) = (input(1)/mass)*((-cos(phi0)*sin(theta0)*sin(psi0))+(sin(phi0)*cos(psi0)));
%
deltaF(8,4) = (input(1)/mass)*((-sin(phi0)*sin(theta0)*sin(psi0))-(cos(phi0)*cos(psi0)));
deltaF(8,5) = (input(1)/mass)*(cos(phi0)*cos(theta0)*sin(psi0));
deltaF(8,6) = (input(1)/mass)*((cos(phi0)*sin(theta0)*cos(psi0))+(sin(phi0)*sin(psi0)));
%
deltaF(9,4) = (input(1)/mass)*(-sin(phi0)*cos(theta0));
deltaF(9,5) = (input(1)/mass)*(-cos(phi0)*sin(theta0));
%
deltaF(10,11) = ((iY-iZ)/iX)*psiDot0;
deltaF(10,12) = ((iY-iZ)/iX)*thetaDot0;
%
deltaF(11,10) = ((iZ-iX)/iY)*psiDot0;
deltaF(11,12) = ((iZ-iX)/iY)*phiDot0;
%
deltaF(12,10) = ((iX-iY)/iZ)*thetaDot0;
deltaF(12,11) = ((iX-iY)/iZ)*phiDot0;
%
phi_k = eye(12) + deltaF*Ts;
% EKF
% Time update
xCap_kMinus = phi_k*(del_x_kp1(:,count));
P_kMinus = phi_k*P_k1*phi_k' + Q;
% Measurement update
K_k = (P_kMinus*H')*inv(H*P_kMinus*H'+R);
%
xCap_k = xCap_kMinus + K_k*(yMeasurement - H*xCap_kMinus);
%
P_k = (eye(12) - K_k*H)*P_kMinus;
%
x_k = xCap_k;
P_k1 = P_k;
%
stateEstimates(:,count) = x_k(1:6);
%
x_k0 = del_x_kp1(:,count);
%
t = Ts:Ts:totT;
%
for k = 1:(length(t)-1)
    %
    x0 = x_k0(1) + e*(-1^count);
    y0 = x_k0(2) + (e*(-1^(count+1)));
    z0 = x_k0(3) + (e*(-1^count));
    xDot0 = x_k0(7) + (e*(-1^(count+1)));
    yDot0 = x_k0(8) + (e*(-1^count));
    zDot0 = x_k0(9) + (e*(-1^(count+1)));
    phi0 = x_k0(4) + (e*(-1^count));
    theta0 = x_k0(5) + (e*(-1^(count+1)));
    psi0 = x_k0(6) + (e*(-1^count));
    phiDot0 = x_k0(10) + (e*(-1^(count+1)));
    thetaDot0 = x_k0(11) + (e*(-1^count));
    psiDot0 = x_k0(12) + (e*(-1^(count+1)));
    %
    simout = sim('EENGquadcopter_control','StartTime',num2str(0),'StopTime',...
        num2str(Ts),'OutputOption','SpecifiedOutputTimes','OutputTimes',...
        num2str(Ts));
    input = simout.Input.Data(end,:)';
    %
    del_x_kp1(:,count+1) = simout.simout.Data(end,:)';
    yMeasurement(:,count+1) = simout.sensorMeasure.Data(end,:)';
    yMeasurement(:,count+1) = yMeasurement(:,count+1) + [sigmaGPS*randn(3,1);...
        sigmaGY*randn(3,1)];
    %
    x = x_k(1);
    y = x_k(2);
    z = x_k(3);
    xDot = x_k(7);
    yDot = x_k(8);
    zDot = x_k(9);
    phi = x_k(4);
    theta = x_k(5);
    psi = x_k(6);
    phiDot = x_k(10);
    thetaDot = x_k(11);
    psiDot = x_k(12);
    % EKF
    deltaF = zeros(12); % Derivative of the non-linear function matrix
    F = zeros(12);
    %
    deltaF(1,7) = 1.0; % xDot
    deltaF(2,8) = 1.0; % yDot
    deltaF(3,9) = 1.0; % zDot
    deltaF(4,10) = 1.0; % phiDot
    deltaF(5,11) = 1.0; % thetaDot
    deltaF(6,12) = 1.0; % psiDot
    %
    deltaF(7,4) = (input(1)/mass)*((-sin(phi)*sin(theta)*cos(psi))+(cos(phi)*sin(psi)));
    deltaF(7,5) = (input(1)/mass)*(cos(phi)*cos(theta)*cos(psi));
    deltaF(7,6) = (input(1)/mass)*((-cos(phi)*sin(theta)*sin(psi))+(sin(phi)*cos(psi)));
    %
    deltaF(8,4) = (input(1)/mass)*((-sin(phi)*sin(theta)*sin(psi))-(cos(phi)*cos(psi)));
    deltaF(8,5) = (input(1)/mass)*(cos(phi)*cos(theta)*sin(psi));
    deltaF(8,6) = (input(1)/mass)*((cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)));
    %
    deltaF(9,4) = (input(1)/mass)*(-sin(phi)*cos(theta));
    deltaF(9,5) = (input(1)/mass)*(-cos(phi)*sin(theta));
    %
    deltaF(10,11) = ((iY-iZ)/iX)*psiDot;
    deltaF(10,12) = ((iY-iZ)/iX)*thetaDot;
    %
    deltaF(11,10) = ((iZ-iX)/iY)*psiDot;
    deltaF(11,12) = ((iZ-iX)/iY)*phiDot;
    %
    deltaF(12,10) = ((iX-iY)/iZ)*thetaDot;
    deltaF(12,11) = ((iX-iY)/iZ)*phiDot;
    %
    phi_k = eye(12) + deltaF*Ts;
    % EKF
    % Time update
    xCap_kMinus = phi_k*(del_x_kp1(:,count+1));
    P_kMinus = phi_k*P_k1*phi_k' + Q;
    % Measurement update
    K_k = (P_kMinus*H')*inv(H*P_kMinus*H'+R);
    %
    xCap_k = xCap_kMinus + K_k*(yMeasurement(:,count+1) - H*xCap_kMinus);
    %
    P_k = (eye(12) - K_k*H)*P_kMinus;
    %
    x_k = xCap_k;
    P_k1 = P_k;
    %
    count = count + 1;
    stateEstimates(:,count) = x_k(1:6);
    %
    x_k0 = del_x_kp1(:,count);
end
%% Plot the results
%
figure
p = plot(t,stateEstimates(3,:),'g',t,yMeasurement(3,:),'b--o',t,del_x_kp1(3,:),'r');
p(1).LineWidth = 2;
title('Z estimate')
xlabel('Time (s)')
ylabel('Z Height (m)')
legend('Estimated State','Noise Measurement','Actual State')
%
figure
p = plot(t,stateEstimates(1,:),'g',t,yMeasurement(1,:),'b--o',t,del_x_kp1(1,:),'r');
p(1).LineWidth = 2;
title('X estimate')
xlabel('Time (s)')
ylabel('X Distance (m)')
legend('Estimated State','Noise Measurement','Actual State')
%
figure
p = plot(t,stateEstimates(2,:),'g',t,yMeasurement(2,:),'b--o',t,del_x_kp1(2,:),'r');
p(1).LineWidth = 2;
title('Y estimate')
xlabel('Time (s)')
ylabel('Y Distance (m)')
legend('Estimated State','Noise Measurement','Actual State')
%
figure
p = plot(t,stateEstimates(4,:),'g',t,yMeasurement(4,:),'b--o',t,del_x_kp1(4,:),'r');
p(1).LineWidth = 2;
title('Phi estimate')
xlabel('Time (s)')
ylabel('Phi angle (rad)')
legend('Estimated State','Noise Measurement','Actual State')
%
figure
p = plot(t,stateEstimates(5,:),'g',t,yMeasurement(5,:),'b--o',t,del_x_kp1(5,:),'r');
p(1).LineWidth = 2;
title('Theta estimate')
xlabel('Time (s)')
ylabel('Theta angle (rad)')
legend('Estimated State','Noise Measurement','Actual State')
%
figure
p = plot(t,stateEstimates(6,:),'g',t,yMeasurement(6,:),'b--o',t,del_x_kp1(6,:),'r');
p(1).LineWidth = 2;
title('Psi estimate')
xlabel('Time (s)')
ylabel('Psi angle (rad)')
legend('Estimated State','Noise Measurement','Actual State')
%