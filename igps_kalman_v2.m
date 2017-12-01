%% Marvelmind iGPS Encoder Kalman Fliter

close all
clear all

%% Import iGPS Data
file_name = 'straight1.log';
file = importdata(file_name,',',2);
file = file.data;

indices = find(file(:,7)==1.85); %Remove beacon locations

file(indices,:) = []; 

t = file(:,1);      %time
x = file(:,5)*1000; 
x = x(1:end-2);     %x
y = file(:,6)*1000;
y = y(1:end-2);     %y
z = file(:,7)*1000;
z = z(1:end-2);     %z (not needed)


%% Variables
X0 = mean(x(1:50)); % X starting point
Y0 = mean(y(1:50)); % Y starting point
Z0 = 0;

theta = 0;
XEnd = mean(x(end-50:end));
YEnd = mean(y(end-50:end));

vLin = 20;       %[mm/s]

fEnc = 50;       %Hz
dtEnc = 1/fEnc;  
encVar = 1; 
SEncR = vLin*dtEnc;
SEncL = vLin*dtEnc;
l = 20;         %wheel base


varBea = 0.1;     
fBea = 50;      %Hz %% Need to take this into account!!!

%% Functions 
% Velocity
X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta)*dtEnc; 
Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta)*dtEnc;
theta_dot_v = @(vr, vl) (1/l)*(vr - vl)*dtEnc;
enc_noise = @(encVar) sqrt(encVar)*randn(1);

% Distances 
X_dot_d = @(dr, dl, theta) (dr/2 + dl/2)*cos(theta); 
Y_dot_d = @(dr, dl, theta) (dr/2 + dl/2)*sin(theta);
theta_dot_d = @(dr, dl) (1/l)*(dr - dl);

%% Initialise 
X = X0;
Y = Y0;

%% Initial Theta
Theta0 = pi - atan((YEnd-Y0)/(X0-XEnd));

%% Simulate and Plot Ideal Data
count = 1;
XIdeal(count,:) = X;
YIdeal(count,:) = Y;
ThetaIdeal(count,:) = Theta0;
%% Generate Simulated Robot Wheel Speeds
while Y<YEnd 
    
    count = count + 1;
    Y = Y + Y_dot_v(vLin, vLin, Theta0);
    X = X + X_dot_v(vLin, vLin, Theta0);
    XIdeal(count,:) = X;
    YIdeal(count,:) = Y;
    ThetaIdeal(count,:) = theta;
    speeds(count,:) = [vLin vLin];
    distances(count,:) = [SEncL SEncR];
    
end

figure
hold
axis equal
plot(XIdeal,YIdeal,'linewidth',2)
scatter(x,y, 'filled')
xlabel('X [mm]')
ylabel('Y [mm]')
title('Ideal Path and iGPS Path','fontsize',17)
legend('Ideal Path','iGPS location','location','southwest')





%% Add Encoder Noise
for i = 1: length(speeds)
encError(i,1) = enc_noise(encVar);
encError(i,2) = enc_noise(encVar);
end

lwv_e = (speeds(:,2) + 1.5*encError(:,2));
rwv_e = (speeds(:,1) + 1.5*encError(:,1));

%% Decimate
factor = round(length(speeds)/length(x),0);

rwv_e = decimate(rwv_e,factor,5);
lwv_e = decimate(lwv_e,factor,5);

%% Reinitialise and Plot

xt = X0;
yt = Y0;
theta = Theta0;


figure
axis equal
for ii = 1:length(rwv_e)
    pose_noise_v(ii,1:4) = [xt, yt, 0, theta];
    
    xt = pose_noise_v(ii,1) + X_dot_v(rwv_e(ii),lwv_e(ii), pose_noise_v(ii,4)); 
    yt = pose_noise_v(ii,2) + Y_dot_v(rwv_e(ii),lwv_e(ii), pose_noise_v(ii,4)); 
    theta = pose_noise_v(ii,4) + (theta_dot_v(rwv_e(ii), lwv_e(ii))/2);
    
end

%% Decimate Encoder Rate to Match GPS
% factor = round(length(speeds)/length(x),0);
% 
% x_dec = decimate(pose_noise_v(:,1),factor,5);
% y_dec = decimate(pose_noise_v(:,2),factor,5);
% theta_dec = decimate(pose_noise_v(:,4),factor,5);
% 
% pose_noise_v = [x_dec,y_dec,theta_dec];
% 
% rwv_e = decimate(rwv_e,factor,5);
% lwv_e = decimate(lwv_e,factor,5);

%%

plot(XIdeal, YIdeal,pose_noise_v(:,1), pose_noise_v(:,2), 'linewidth',2)
xlabel('X [mm]')
ylabel('Y [mm]')
%axis equal
title('Ideal Path and Encoder Dead Reckoning','fontsize',17)
legend('Ideal Path','Dead Reckoning','location','southwest')

%% EKF

R = [varBea 0 0 
     0 varBea 0 
     0 0 varBea ]; 
    
V = eye(3);

Q = [encVar 0
    0   encVar];

P_ini = ones(4);
P = P_ini;
X_kal = zeros(length(x),4);

% Initial Guess of State

X_hat_prev = [X0, Y0, Z0, Theta0]; 


X_hat_minus_plot = zeros(length(x),4);



for i = 1:length(x)

       
    X_hat_minus = [X_hat_prev(1) + (1/2)*(rwv_e(i)+lwv_e(i))*cos(X_hat_prev(4))*dtEnc
                   X_hat_prev(2) + (1/2)*(rwv_e(i)+lwv_e(i))*sin(X_hat_prev(4))*dtEnc
                   X_hat_prev(3)
                   X_hat_prev(4) + (lwv_e(i)-rwv_e(i))*dtEnc/l   ];
               
    X_hat_prev = X_hat_minus; 
               
    X_hat_minus_log(i,:) = X_hat_minus;
end;
figure;
subplot(3,1,1);
plot(pose_noise_v(:,1), pose_noise_v(:,2), 'linewidth',2)
hold on;
subplot(3,1,2);
plot(X_hat_minus_log(:,1),X_hat_minus_log(:,2));
subplot(3,1,3);
plot(XIdeal, YIdeal);
xlabel('X [mm]')
ylabel('Y [mm]')
legend('x_hat', 'Ideal', 'Dead Reckoning');

%% Reinitialise

X_hat_prev = [X0, Y0, Z0, Theta0]; 


X_hat_minus_plot = zeros(length(x),4);


for i = 1:length(x)

       
    X_hat_minus = [X_hat_prev(1) + (1/2)*(rwv_e(i)+lwv_e(i))*cos(X_hat_prev(4))
                   X_hat_prev(2) + (1/2)*(rwv_e(i)+lwv_e(i))*sin(X_hat_prev(4))
                   X_hat_prev(3)
                   X_hat_prev(4) + (lwv_e(i)-rwv_e(i))*dtEnc/l   ];
            
    % Process Jacobian
    
    A = [1  0   0   -(1/2)*(rwv_e(i) + lwv_e(i))*dtEnc*sin(X_hat_prev(4))
        0   1   0  (1 /2)*(rwv_e(i) + lwv_e(i))*dtEnc*cos(X_hat_prev(4))
        0   0   1   0
        0   0   0   1];
    
    W = [(1/2)*cos(X_hat_prev(4))*dtEnc  (1/2)*sin(X_hat_prev(4))*dtEnc
        (1/2)*cos(X_hat_prev(4))*dtEnc  (1/2)*sin(X_hat_prev(4))*dtEnc
        0   0
        (dtEnc/l)     (dtEnc/l)];
      
    
    % Measurment Jacobian
    H = [   1 0 0 0
            0 1 0 0
            0 0 1 0    ];
    
        
    z = [x(i)
         y(i)
         0];
  
    
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


figure
plot(X_kal(:,1),X_kal(:,2))



figure
hold on

plot(pose_noise_v(:,1), pose_noise_v(:,2),'LineWidth',2)
plot(X_kal(:,1),X_kal(:,2),'LineWidth',2)
plot(XIdeal, YIdeal,'LineWidth',2)
scatter(x,y, 'filled')
title('Kalman Filter Positioning Comparasion')
xlabel('X Position')
ylabel('Y Position')
legend('Dead Reckoning','Kalman Filter','Straight Line','Location','Bestoutside', 'Location', 'northeast')


