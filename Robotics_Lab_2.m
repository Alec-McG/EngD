close all
clear all
%% Variables
vMax = 20; %mm/s
x_0 = -300;
y_0 = -300;
theta_0 = 0;
l=100;
r = 10;
x = x_0;    %initialise
y = y_0;
theta = theta_0; 
t = 0;
phiDot_l = vMax;
phiDot_r = vMax;
fEnc = 10; %Hz
dtEnc = 1/fEnc; %Encoder Timestep
encVar = 1; %mm/s2
%% beacons
b1 = [2500,2500,1000];
b2 = [-2500, 2500,1000];
b3 = [2500,-2500,1000];
varBea = 23;
fBea = 5; %Hz
%% Functions
x_dot = @(vr, vl) (vr./2 + vl./2)*dtEnc; 
y_dot = @(vr, vl) (vr./2 + vl./2)*dtEnc;
theta_dot = @(vr, vl) (r/l)*(vr - vl)*dtEnc; 
%% Get Wheel Speeds and Plot Ideal Location
count = 0;

for ii = 1:7 
    
    x_0 = x;
    theta_0 = 0;
    y_0 = y;
    
    if ii < 7
        if mod(ii,2)==1
            
            while abs(x-x_0)<2000
            count = count + 1;
            x = x + x_dot(phiDot_l, phiDot_r);
            x_pos(count,:) = x;
            wheel_speeds(count,:) = [phiDot_l phiDot_r];
            y_pos(count,:) = y;
            theta_pos(count,:) = theta;
        
            end

            while theta - theta_0 < pi/2
            count = count + 1;
            theta = theta + theta_dot(vMax, -vMax);
            theta_pos(count,:) = theta;
            wheel_speeds(count,:) = [vMax -vMax];
            x_pos(count,:) = x;
            y_pos(count,:) = y;
            end
        
            
            while (y-y_0)<300
            count = count + 1;
            y = y + y_dot(vMax, vMax);
            y_pos(count,:) = y;
            wheel_speeds(count,:) = [vMax vMax];
            x_pos(count,:) = x;
            theta_pos(count,:) = theta;
            end
            
            while theta - theta_0 < pi
            count = count + 1;
            theta = theta + theta_dot(vMax, -vMax);
            theta_pos(count,:) = theta;
            wheel_speeds(count,:) = [vMax -vMax];
            x_pos(count,:) = x;
            y_pos(count,:) = y;
            end
        
        else
       
            while abs(x-x_0)<2000
            count = count + 1;
            x = x - x_dot(vMax, vMax);
            x_pos(count,:) = x;
            wheel_speeds(count,:) = [vMax vMax];
            y_pos(count,:) = y;
            theta_pos(count,:) = theta;
            end

            while theta - theta_0 > pi/2
            count = count + 1;
            theta = theta + theta_dot(-vMax, vMax);
            theta_pos(count,:) = theta;
            wheel_speeds(count,:) = [-vMax vMax];
            y_pos(count,:) = y;
            x_pos(count,:) = x;
            end
       
            
            
            while abs(y-y_0)<300
            count = count + 1;
            y = y + y_dot(vMax, vMax);
            y_pos(count,:) = y;
            wheel_speeds(count,:) = [vMax vMax];
            x_pos(count,:) = x;
            theta_pos(count,:) = theta;
            end
            
            while theta - theta_0 > 0
            count = count + 1;
            theta = theta + theta_dot(-vMax, vMax);
            theta_pos(count,:) = theta;
            wheel_speeds(count,:) = [-vMax vMax];
            y_pos(count,:) = y;
            x_pos(count,:) = x;
            end
        end
    if ii == 7
        
    end
    else 
        while abs(x-x_0)<2000
            count = count + 1;
            x = x + x_dot(phiDot_l, phiDot_r);
            x_pos(count,:) = x;
            wheel_speeds(count,:) = [phiDot_l phiDot_r];
            y_pos(count,:) = y;
            theta_pos(count,:) = theta;
        end
    end
end

view(3)
figure(1);
hold on
plot3(x_pos, y_pos, zeros(length(x_pos),1));
plot3(b1(1), b1(2), b1(3), '*r');
plot3(b2(1), b2(2), b2(3), '*r')
plot3(b3(1), b3(2), b3(3), '*r')
hold off

% figure(2)
% subplot(1,2,1);
% plot(wheel_speeds(:,1));
% subplot(1,2,2)
% plot(wheel_speeds(:,2));

%% Add Add Wheel Encoder Noise (Dead-Reckoning)`
enc_noise = @(encVar) sqrt(encVar)*randn(1);
x_dot = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta)*dtEnc; 
y_dot = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta)*dtEnc;
theta_dot = @(vr, vl) (r/l)*(vr - vl)*dtEnc;

%% get speed errors
for i = 1: length(wheel_speeds)
encError(i,1) = enc_noise(encVar);
end
for i = 1: length(wheel_speeds)
encError(i,2) = enc_noise(encVar);
end
%% New wheel Speeds
lw_e = (wheel_speeds(:,2) + encError(:,2));
rw_e = (wheel_speeds(:,1) + encError(:,1));


%% reinitialise
xt = -300;
yt = -300;
theta = 0;
%pose = zeros(length(lw_e),3);

figure(3)
for ii = 1:length(lw_e)
    pose(ii,1:4) = [xt, yt, 0, theta];
    
    xt = pose(ii,1) + x_dot(rw_e(ii),lw_e(ii), pose(ii,4)); 
    yt = pose(ii,2) + y_dot(rw_e(ii),lw_e(ii), pose(ii,4)); 
    theta = pose(ii,4) + theta_dot(rw_e(ii), lw_e(ii));
        
end

plot(pose(:,1), pose(:,2))




%% EKF

%% From lectures

abs_pos = [x_pos, y_pos, zeros(length(x_pos),1)];

h_func = @(x, y, b, varBea) sqrt((x - b(1)).^2 +(y - b(2)).^2 + (0 - b(3)).^2) + sqrt(varBea)*randn(1);
    
R = [varBea 0 0 
     0 varBea 0 
     0 0 varBea ]; 
    
V = eye(4);

Q = [encVar 0
    0   encVar];

P_ini = diag([2*encVar, 2*encVar, 0,0]);

%P_initial = R;
P = P_ini;
X_kal = zeros(length(x_pos),3);

X_pre = pose(1,1:4)';

figure
for i = 1:length(x_pos)
   

    % Process Jacobian
    
    A = [1  0   0   -0.5*(wheel_speeds(i,1) + wheel_speeds(i,2))*sin(pose(i,3))
        0   1   0  0.5*(wheel_speeds(i,1) + wheel_speeds(i,2))*cos(pose(i,3)) 
        0   0   1   0
        0   0   0   1];
    
    W = [0.5*cos(pose(i,3))  0.5*sin(pose(i,3))
        0.5*cos(pose(i,3))  0.5*sin(pose(i,3))
        0   0
        1/l     1/l];
      
    z = [h_func(x_pos(i), y_pos(i), b1, varBea) h_func(x_pos(i), y_pos(i), b2, varBea) h_func(x_pos(i), y_pos(i), b3, varBea)   0]';
    %h = [h_func(x_pos(i), y_pos(i), b1, 0) h_func(x_pos(i), y_pos(i), b2, 0) h_func(x_pos(i), y_pos(i), b3, 0)]';
    
    L1 = sqrt((x_pos(i)-b1(1))^2 +(y_pos(i)-b1(2))^2 + (0 - b1(3))^2) ;
    L2 = sqrt((x_pos(i)-b2(1))^2 +(y_pos(i)-b2(2))^2 + (0 - b2(3))^2);
    L3 = sqrt((x_pos(i)-b3(1))^2 +(y_pos(i)-b3(2))^2 + (0 - b3(3))^2);
    
    % Measurment Jacobian
    H = [   (x_pos(i)-b1(1))/ L1    (y_pos(i) -b1(2))/L1   (0 - b1(3))/L1   0
            (x_pos(i)-b2(1))/ L2    (y_pos(i) -b2(2))/L2   (0 - b2(3))/L2   0
            (x_pos(i)-b3(1))/ L3    (y_pos(i) -b3(2))/L3   (0 - b3(3))/L3   0
                0   0   0   0];
    
   % Kalman Filter Implimentation
   % Kalman Gain
   
   P = A*P*A' + P_ini;
   K = P.*H'.*inv(H.*P.*H' +V.*P_ini.*V');
   
   
   
   
   %EKF Estiamte
   X_pre = X_pre + K*(z - pose(i,1:4)'); % h? meaurment and model???    
   
   X_kal(i,1:4) = X_pre';
   
   P = (eye(4)-K*H)*P; 
  
   plot(X_kal(:,1),X_kal(:,2))

    
end


%% Measurment Jacobians




