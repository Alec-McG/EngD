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

X_Previous = X0;
Y_Previous = Y0;
Z_Previous = Z0;

theta = pi/2;
XEnd = mean(x(end-50:end));
YEnd = mean(y(end-50:end));

vLin = 20;       %[mm/s]

fEnc = 10;       %Hz
dtEnc = 1/fEnc;  
encVar = 0.1; 

l = 20;         %wheel base

fGPS = 1;       %Hz
varBea = 20;     

r = 10;     %Wheel Radius
C = 2*pi*r; %Wheel circumfrence

%% Functions 
% Velocity

theta_dot_v = @(vr, vl) (1/l)*(vr - vl)*dtEnc;
X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta + (theta_dot_v(vr, vl)/2))*dtEnc; 
Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta + (theta_dot_v(vr, vl)/2))*dtEnc;
% 
% X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta)*dtEnc; 
% Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta)*dtEnc;

enc_noise = @(encVar) sqrt(encVar)*randn(1);

% Distances 
% X_dot_d = @(dr, dl, theta) (dr/2 + dl/2)*cos(theta); 
% Y_dot_d = @(dr, dl, theta) (dr/2 + dl/2)*sin(theta);
% theta_dot_d = @(dr, dl) (1/l)*(dr - dl);

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

lwv_e = (speeds(:,2) + encError(:,2));
rwv_e = (speeds(:,1) + encError(:,1));

%% Simulate GPS data
x_sim_gps = zeros(1,length(XIdeal));
y_sim_gps = zeros(1,length(XIdeal));

tt = linspace(0,2*pi,length(XIdeal));

for q = 1:length(XIdeal)
    
    x_sim_gps(q) = XIdeal(q) + randn()*(varBea) + 0*sin(tt(q)/2); 
    y_sim_gps(q) = YIdeal(q) + randn()*(varBea);
    
end

factor = fEnc/fGPS;
x_sim_gps = decimate(x_sim_gps,factor,5);
y_sim_gps = decimate(y_sim_gps,factor,5);

%% Reinitialise and Plot

xt = X0;
yt = Y0;
theta = Theta0;

figure
axis equal
for ii = 1:length(speeds)
    pose_noise_v(ii,1:4) = [xt, yt, 0, theta];
    
    xt = pose_noise_v(ii,1) + X_dot_v(rwv_e(ii),lwv_e(ii), pose_noise_v(ii,4)); 
    yt = pose_noise_v(ii,2) + Y_dot_v(rwv_e(ii),lwv_e(ii), pose_noise_v(ii,4)); 
    theta = pose_noise_v(ii,4) + (theta_dot_v(rwv_e(ii), lwv_e(ii)));
    
end


%% Decimate and Convert Encoders and Units 

factor = fEnc/fGPS;
lwv_e = decimate(lwv_e,factor,5)/dtEnc;
rwv_e = decimate(rwv_e,factor,5)/dtEnc;

%% Simulated GPS

figure();
scatter(x,y);

hold on;
scatter(x_sim_gps,y_sim_gps)
title('Example gps data for x-coordinate');

%%
figure;
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

P_ini = zeros(4,4);
P = P_ini;
X_kal = zeros(length(x_sim_gps),4);

% Initial Guess of State

X_hat_prev = [X0, Y0, Z0, Theta0]; 


X_hat_minus_plot = zeros(length(x),4);

X_hat_prev(1) = X0;
X_hat_prev(2) = Y0;
X_hat_prev(3) = 0;
X_hat_prev(4) = Theta0;

X_kal(1,1:4) = [X_hat_prev(1) X_hat_prev(2) X_hat_prev(3) X_hat_prev(4)];

count1 = 0;
count2 = 0;

for i = 1:length(rwv_e)
          
    X_hat_minus = [X_hat_prev(1) + X_dot_v(rwv_e(i),lwv_e(i),X_hat_prev(4))
                   X_hat_prev(2) + Y_dot_v(rwv_e(i),lwv_e(i),X_hat_prev(4))
                   X_hat_prev(3)
                   X_hat_prev(4) + theta_dot_v(rwv_e(i),lwv_e(i))   ];
    
%     %If the GPS coordinate is the same as the previous then do not use a
%     %correction factor, rely on the dead reckoning only
%      if((x_sim_gps(i) == X_Previous)  && (y_sim_gps(i) == Y_Previous))
%     %Only DeadReconing
%         X_hat = X_hat_minus;
%         X_hat_prev = X_hat;
%         X_kal(i,1:4) = X_hat(1:4)';
%   
%     else

    
 count1 = count1 + 1;
 
 % Use previous GPS coordinates for correction
 
%     if((x_sim_gps(i) == X_Previous)  && (y_sim_gps(i) == Y_Previous))
%         
%     count2 = count2 + 1;
% 
%     
%     A = [1  0   0   -(1/2)*(rwv_e(i) + lwv_e(i))*sin(X_hat_prev(4))
%         0   1   0  (1/2)*(rwv_e(i) + lwv_e(i))*cos(X_hat_prev(4))
%         0   0   1   0
%         0   0   0   1];
%     
%     W = [(1/2)*cos(X_hat_prev(4))  (1/2)*cos(X_hat_prev(4))
%         (1/2)*sin(X_hat_prev(4)) (1/2)*sin(X_hat_prev(4))
%         0   0
%         ( 1/l)     (1/l)];
%       
%     
%     % Measurment Jacobian
%     H = [   1 0 0 0
%             0 1 0 0
%             0 0 1 0    ];
%     
%         
%     z = [x_sim_gps(i-(count1-count2)) % Use the previous X position
%          X_hat_prev(2)        % Use dead reckoning for y coordinate propogation
%          0];
%   
%     
%     h_0 = H*X_hat_minus; %% from online presentation???
%     
%    % Kalman Filter Implimentation
%      
%    P = A*P*A' + W*Q*W'; 
%   
%    K = P*H'*inv(H*P*H' + V*R*V'); %%%%%
%    
%    %EKF Correction
%    X_hat = X_hat_minus  + K*(z - h_0);
%    
%    X_kal(i,1:4) = X_hat(1:4)';
%    
%    P = (eye(4)-K*H)*P; 
%    
%    X_hat_prev = X_hat;
% 
%     
% else 
    

    
    %Perform Kalman filter
    % Process Jacobian
    
    A = [1  0   0   -(dtEnc/2)*(rwv_e(i) + lwv_e(i))*sin(X_hat_prev(4))
        0   1   0  (dtEnc/2)*(rwv_e(i) + lwv_e(i))*cos(X_hat_prev(4))
        0   0   1   0
        0   0   0   1];
    
    W = [(dtEnc/2)*cos(X_hat_prev(4))  (dtEnc/2)*cos(X_hat_prev(4))
        (dtEnc/2)*sin(X_hat_prev(4)) (dtEnc/2)*sin(X_hat_prev(4))
        0   0
        ( -dtEnc/l)     (dtEnc/l)];
      
      
    
    % Measurment Jacobian
    H = [   1 0 0 0
            0 1 0 0
            0 0 1 0    ];
    
        
    z = [x_sim_gps(i)
         y_sim_gps(i)
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
   
%     end;
%     
   X_Previous = x(i);
   Y_Previous = y(i);
   %Z_Previous = z(i);
    
end

% 
% figure
% plot(X_kal(:,1),X_kal(:,2))



figure
hold on

plot(pose_noise_v(:,1), pose_noise_v(:,2),'r')
scatter(X_kal(:,1),X_kal(:,2),'g','filled')
plot(X_kal(:,1),X_kal(:,2),'g')
plot(XIdeal, YIdeal,'b')



% plot(X_kal(1:10:length(X_kal),1),X_kal(1:10:length(X_kal),2),'m')
scatter(x_sim_gps,y_sim_gps, 'filled')
title('Kalman Filter Positioning Comparasion','FontSize', 20)
xlabel('X Position')
ylabel('Y Position')
legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','GPS Measurment','Bestoutside', 'Location', 'northeast')

grid on

% fit = polyfit(X_kal(:,1),X_kal(:,2),1);
% plot(polyval(fit,X_kal(:,1)))
