%% Marvelmind iGPS Encoder Kalman Fliter

close all
clear all

%% Variables
l = 50;                     % Wheel base
r = 10;                     % Wheel radius
C = 2*pi*r;                 % Wheel circumference
v = 100;                    % 20mm/s max linear velocity
w = 2*v/l;                  % Angular speed when rotating
fs_IMU = 5 ;                % 5 Hz  
fs_enc = 10;                % 10 Hz  
dtEnc = 1/fs_enc;           % Encoder timestep
encVar = 0.5;               % Encoder Variance
sim_length = 10000;         % Straight line length
t_len = sim_length/v;       % Time moving 

gyroThetaVar = 2*pi/180;    % Variance of Gyroscope measurment

X0 = 0;
Y0 = 0;
Theta0 = 90*pi/180;


%% Functions 
% Velocity
theta_dot_v = @(vr, vl) (1/l)*(vr - vl)*dtEnc;
X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta + (theta_dot_v(vr, vl)/2))*dtEnc; 
Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta + (theta_dot_v(vr, vl)/2))*dtEnc;

%X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta)*dtEnc; 
%Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta)*dtEnc;

enc_noise = @(encVar) sqrt(encVar)*randn(1);


%% Initialise 
X = X0;
Y = Y0;
Theta = Theta0;


%% Simulate and Plot Ideal Data
count = 1;
XIdeal(count,:) = X;
YIdeal(count,:) = Y;
ThetaIdeal(count,:) = Theta0;
%% Generate Simulated Robot Wheel Speeds

while Y < sim_length
    
    count = count + 1;
    Y = Y + Y_dot_v(v, v, Theta0);
    X = X + X_dot_v(v, v, Theta0);
    XIdeal(count,:) = X;
    YIdeal(count,:) = Y;
    ThetaIdeal(count,:) = Theta0;
    speeds(count,:) = [v v];
        
end

figure(1)
plot(XIdeal,YIdeal)
%axis equal
xlabel('X [mm]')
ylabel('Y [mm]')
title('Ideal Path','fontsize',17)
legend('Ideal Path','location','southwest')




for ii=1:10
%% Dead Reckoning

for j = 1: length(speeds)
    encError(j,1) = enc_noise(encVar);
    encError(j,2) = enc_noise(encVar);
end

lwv_e = (speeds(:,2) + encError(:,2));
rwv_e = (speeds(:,1) + encError(:,1));

% Reinitialise and Plot
xt = X0;
yt = Y0;
theta = Theta0;

figure
%axis equal
%% Simulate Measurements and Encoder Data

for jj = 1:length(rwv_e)
    pose_noise_v(jj,1:4) = [xt, yt, 0, theta];
    
    xt = pose_noise_v(jj,1) + X_dot_v(rwv_e(jj),lwv_e(jj), pose_noise_v(jj,4)); 
    yt = pose_noise_v(jj,2) + Y_dot_v(rwv_e(jj),lwv_e(jj), pose_noise_v(jj,4)); 
    theta = pose_noise_v(jj,4) + (theta_dot_v(rwv_e(jj), lwv_e(jj))/2);
    thetaGyroMeas(jj) = Theta0 + randn(1)*gyroThetaVar; 
end

% plot(XIdeal, YIdeal,pose_noise_v(:,1), pose_noise_v(:,2), 'linewidth',2)
% xlabel('X [mm]')
% ylabel('Y [mm]')
% %axis equal
% title('Ideal Path and Encoder Dead Reckoning','fontsize',17)
% legend('Ideal Path','Dead Reckoning','location','southwest')


thetaGyroMeas = thetaGyroMeas';

%% EKF

R = [gyroThetaVar^2]; 
    
V = eye(1);

Q = [encVar^2 0
    0   encVar^2];
  
P_ini = zeros(4,4);
P = P_ini;
X_kal = zeros(length(pose_noise_v),4);

% Initial Guess of State

X_hat_prev = [X0, Y0, 0, Theta0]; 

X_hat_minus_plot = zeros(length(pose_noise_v),4);

X_hat_prev(1) = X0;
X_hat_prev(2) = Y0;
X_hat_prev(3) = 0;
X_hat_prev(4) = Theta0;

X_kal = zeros(length(rwv_e),4);

X_kal(1,1:4) = [X_hat_prev(1) X_hat_prev(2) X_hat_prev(3) X_hat_prev(4)];  
    
for i = 1:length(rwv_e)
      
       
    X_hat_minus = [X_hat_prev(1) + X_dot_v(rwv_e(i),lwv_e(i),X_hat_prev(4))
                   X_hat_prev(2) + Y_dot_v(rwv_e(i),lwv_e(i),X_hat_prev(4))
                   X_hat_prev(3)
                   X_hat_prev(4) + theta_dot_v(rwv_e(i),lwv_e(i))   ];
    

    
    A = [1  0   0   -(dtEnc/2)*(rwv_e(i) + lwv_e(i))*sin(X_hat_prev(4))
        0   1   0  (dtEnc/2)*(rwv_e(i) + lwv_e(i))*cos(X_hat_prev(4))
        0   0   1   0
        0   0   0   1];
    
    W = [(dtEnc/2)*cos(X_hat_prev(4))  (dtEnc/2)*cos(X_hat_prev(4))
        (dtEnc/2)*sin(X_hat_prev(4)) (dtEnc/2)*sin(X_hat_prev(4))
        0   0
        ( dtEnc/l)     (dtEnc/l)];
      
    
    % Measurment Jacobian
    H = [   0 0 0 1   ];
    
        
    z = [thetaGyroMeas(i)]; % Use the previous X position
         

  
    
    h_0 = H*X_hat_minus; %% from online presentation???
    
   % Kalman Filter Implimentation
     
   P = A*P*A' + W*Q*W'; 
  
   K = P*H'*inv(H*P*H' + V*R*V'); %%%%%
   
   %EKF Correction
   X_hat = X_hat_minus  + K*(z - h_0);
   
   X_kal(i,1:4) = X_hat(1:4)';
   
   P = (eye(4)-K*H)*P; 
   
   X_hat_prev = X_hat;

    
end

figure(ii+1)
hold on
plot(pose_noise_v(:,1), pose_noise_v(:,2),'r')
scatter(X_kal(:,1),X_kal(:,2),'g')
plot(X_kal(:,1),X_kal(:,2),'g')
plot(XIdeal, YIdeal,'b')

title('Kalman Filter Positioning Comparasion','FontSize', 20)
xlabel('X Position [mm]')
ylabel('Y Position [mm]')
legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','GPS Measurment','Bestoutside', 'Location', 'northeast')
hold off

grid on
pause(0.1)
end