%% Unscented Kalman Filter
%
Ts = 0.025;
totT = 2.5;
sigmaGPS = 0.25;
sigmaGY = 0.1;
e = 1e-6;
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
count = 1;
%
t = Ts:Ts:totT;
%
P_k1 = eye(12);
%
Q = diag(e*ones(1,12));
%
R = eye(6);
R(1,1) = sigmaGPS*R(1,1);
R(2,2) = sigmaGPS*R(2,2);
R(3,3) = sigmaGPS*R(3,3);
R(4,4) = sigmaGY*R(4,4);
R(5,5) = sigmaGY*R(5,5);
R(6,6) = sigmaGY*R(6,6);
%
stateEstimates = zeros(6,24);
%
x_k0 = [x0;y0;z0;phi0;theta0;psi0;...
    xDot0;yDot0;zDot0;phiDot0;thetaDot0;psiDot0];
%
n = length(x_k0);
%
for k = 1:length(t)
    %
    x0 = x_k0(1)+e*(-1^(count+1));
    y0 = x_k0(2)+e*(-1^count);
    z0 = x_k0(3)+e*(-1^(count+1));
    xDot0 = x_k0(7)+e*(-1^count);
    yDot0 = x_k0(8)+e*(-1^(count+1));
    zDot0 = x_k0(9)+e*(-1^count);
    phi0 = x_k0(4)+e*(-1^(count+1));
    theta0 = x_k0(5)+e*(-1^count);
    psi0 = x_k0(6)+e*(-1^(count+1));
    phiDot0 = x_k0(10)+e*(-1^count);
    thetaDot0 = x_k0(11)+e*(-1^(count+1));
    psiDot0 = x_k0(12)+e*(-1^count);
    %
    simout = sim('EENGquadcopter_control','StartTime',num2str(0),'StopTime',...
            num2str(Ts),'OutputOption','SpecifiedOutputTimes','OutputTimes',...
            num2str(Ts));
    %
    input = simout.Input.Data(end,:)';
    %
    del_x_kp1(:,count) = simout.simout.Data(end,:)';
    yMeasurement(:,count) = simout.sensorMeasure.Data(end,:)';
    yMeasurement(:,count) = yMeasurement(:,count) + [sigmaGPS*randn(3,1);...
        sigmaGY*randn(3,1)];
    % Time update
    M = chol(P_k1,'upper');
    %
    xBar_i_p = sqrt(n)*[M -M];
    %
    xI_p = x_k0 + xBar_i_p;
    %
    for i = 1:(2*n)
        %
        x0 = xI_p(1,i);
        y0 = xI_p(2,i);
        z0 = xI_p(3,i);
        phi0 = xI_p(4,i);
        theta0 = xI_p(5,i);
        psi0 = xI_p(6,i);
        xDot0 = xI_p(7,i);
        yDot0 = xI_p(8,i);
        zDot0 = xI_p(9,i);
        phiDot0 = xI_p(10,i);
        thetaDot0 = xI_p(11,i);
        psiDot0 = xI_p(12,i);
        %
        simout_Temp = sim('EENGquadcopter_control','StartTime',num2str(0),'StopTime',...
            num2str(Ts),'OutputOption','SpecifiedOutputTimes','OutputTimes',...
            num2str(Ts));
        %
        x_ki(:,i) = simout_Temp.simout.Data(end,:)';
    end
    %
    xCap_Minus = (1/(2*n))*sum(x_ki,2);
    %
    PMinus = zeros(12);
    %
    for i = 1:(2*n)
       PMinus = PMinus + (x_ki(:,i)-xCap_Minus)*(x_ki(:,i)-xCap_Minus)' + Q; 
    end
    %
    PMinus = (1/(2*n)).*PMinus;
    %
    % Measurement update
    M = chol(PMinus,'upper');
    %
    xBar_i = sqrt(n)*[M -M];
    %
    xI_p = xCap_Minus + xBar_i;
    %
    for i = 1:(2*n)
        %
        x0 = xI_p(1,i);
        y0 = xI_p(2,i);
        z0 = xI_p(3,i);
        phi0 = xI_p(4,i);
        theta0 = xI_p(5,i);
        psi0 = xI_p(6,i);
        xDot0 = xI_p(7,i);
        yDot0 = xI_p(8,i);
        zDot0 = xI_p(9,i);
        phiDot0 = xI_p(10,i);
        thetaDot0 = xI_p(11,i);
        psiDot0 = xI_p(12,i);
        %
        simout_Temp = sim('EENGquadcopter_control','StartTime',num2str(0),'StopTime',...
            num2str(Ts),'OutputOption','SpecifiedOutputTimes','OutputTimes',...
            num2str(Ts));
        %
        yMeasure(:,i) = simout_Temp.sensorMeasure.Data(end,:)';
    end
    %
    yBar = (1/(2*n))*sum(yMeasure,2);
    %
    Pxy = zeros(12,6);
    Py = zeros(6);
    %
    for i = 1:(2*n)
       Pxy = Pxy + (xI_p(:,i)-xCap_Minus)*(yMeasure(:,i)-yBar)';
       Py = Py + (yMeasure(:,i)-yBar)*(yMeasure(:,i)-yBar)' + R;
    end
    %
    Pxy = (1/(2*n)).*Pxy;
    Py = (1/(2*n)).*Py;
    %
    K = Pxy*inv(Py);
    %
    xCap = xCap_Minus + K*(yMeasurement(:,count) - yBar);
    %
    P_k = PMinus - K*Pxy';
    %
    stateEstimates(:,count) = xCap(1:6);
    %
    x_k0 = del_x_kp1(:,count);
    %
    count = count + 1;
end
%
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