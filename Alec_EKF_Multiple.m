close all
clear all
%% Notes
%   In each option IMU and GPS sampling rates are the same and encoder
%   simulation is decimated to match.
%
%   Encoder + IMU + GPS simulation does not take into account different
%   sampling rates. The lowest sampling rate sbhould be used for decimation 
%
%% Pick An Option Below

Option = 1;

%   1 = 2 Encoder + GPS
%   2 = 2 Encoder + IMU
%   3 = 2 Encoder + IMU + GPS
%   4 = 1 Encoder + IMU 


%% Variables
l = 50;                     % Wheel base
r = 10;                     % Wheel radius
C = 2*pi*r;                 % Wheel circumference
v = 100;                    % 20mm/s max linear velocity
w = 2*v/l;                  % Angular speed when rotating
fIMU = 1;                   % 1 Hz  
fGPS = 1;                   % 1 Hz 
fEnc = 10;                  % 5 Hz  

dtEnc = 1/fEnc;             % Encoder timestep

sim_length = 10000;         % Straight line length
t_len = sim_length/v;       % Time moving 

encStd = 0.1;               % Standard Deviations
GPSStd = 200;  
IMUStd = 0.1*pi/180;          

sinamp = 0;                %Amplitude of GPS X sine wave

X0 = 0;
Y0 = 0;
Z0 = 0;
Theta0 = 90*pi/180;

IMU_points = t_len/(fIMU);
gps_points = t_len/(fGPS);

%% Functions 
% Velocity
theta_dot_v = @(vr, vl) (1/l)*(vr - vl)*dtEnc;
X_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*cos(theta + (theta_dot_v(vr, vl)/2))*dtEnc; 
Y_dot_v = @(vr, vl, theta) (vr/2 + vl/2)*sin(theta + (theta_dot_v(vr, vl)/2))*dtEnc;
enc_noise = @(encVar) (encVar)*randn(1);

%% Initialise 
X = X0;
Y = Y0;
Theta = Theta0;

%% Simulate and Plot Ideal Data
count = 1;
XIdeal(count,:) = X;
YIdeal(count,:) = Y;
ThetaIdeal(count,:) = Theta0;
speeds(count,:)= [v v];

while Y < sim_length
    
    count = count + 1;
    Y = Y + Y_dot_v(v, v, Theta0);
    X = X + X_dot_v(v, v, Theta0);
    XIdeal(count,:) = X;
    YIdeal(count,:) = Y;
    ThetaIdeal(count,:) = Theta0;
    speeds(count,:) = [v v];
        
end


%% Preset Sizes
encError = zeros(length(speeds),2);
pose_noise = zeros(length(speeds),4);
thetaGyroMeas = zeros(IMU_points,1);
x_sim_gps = zeros(1,length(gps_points));
y_sim_gps = zeros(1,length(gps_points));
X_kal = zeros(IMU_points,4);

for pp =1:4
            %% Add Encoder Noise
        for i = 1: length(speeds)
            encError(i,1) = enc_noise(encStd);
            encError(i,2) = enc_noise(encStd);
        end

        lwv_e = (speeds(:,2) + encError(:,2));
        rwv_e = (speeds(:,1) + encError(:,1));

        
        %% Reinitialise and Plot

        xt = X0;
        yt = Y0;
        theta = Theta0;

        for ii = 1:length(speeds)
            pose_noise(ii,1:4) = [xt, yt, 0, theta];

            xt = pose_noise(ii,1) + X_dot_v(rwv_e(ii),lwv_e(ii), pose_noise(ii,4)); 
            yt = pose_noise(ii,2) + Y_dot_v(rwv_e(ii),lwv_e(ii), pose_noise(ii,4)); 
            theta = pose_noise(ii,4) + (theta_dot_v(rwv_e(ii), lwv_e(ii)));

        end
        
        if Option == 1 %Encoder + GPS
            %% Simulate GPS data
            x_sim_gps = zeros(1,length(XIdeal));
            y_sim_gps = zeros(1,length(XIdeal));

            tt = linspace(0,2*pi,length(XIdeal));

            for q = 1:length(XIdeal)

                x_sim_gps(q) = XIdeal(q) + randn()*(GPSStd) + sinamp*sin(tt(q)/2); 
                y_sim_gps(q) = YIdeal(q) + randn()*(GPSStd);

            end

            factor = fEnc/fGPS;
            x_sim_gps = decimate(x_sim_gps,factor,5);
            y_sim_gps = decimate(y_sim_gps,factor,5);

            %% Decimate and Convert Encoders and Units 

            factor = fEnc/fGPS;
            lwv_e = decimate(lwv_e,factor,5)/dtEnc;
            rwv_e = decimate(rwv_e,factor,5)/dtEnc;

            %% EKF

            R = [GPSStd 0 0 
                0 GPSStd 0 
                0 0 GPSStd ]; 

            V = eye(3);

            Q = [encStd 0
                0   encStd];

            P_ini = zeros(4,4);
            P = P_ini;
            X_kal = zeros(length(x_sim_gps),4);

            % Initial Guess of State

            X_hat_prev = [X0, Y0, Z0, Theta0]; 

            X_hat_prev(1) = X0;
            X_hat_prev(2) = Y0;
            X_hat_prev(3) = 0;
            X_hat_prev(4) = Theta0;

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

                K = P*H'*(H*P*H' + V*R*V')^-1; %%%%%

                %EKF Correction
                X_hat = X_hat_minus  + K*(z - h_0);

                X_kal(i+1,1:4) = X_hat(1:4)';

                P = (eye(4)-K*H)*P; 

                X_hat_prev = X_hat;

end 
        
    elseif Option == 2 % Encoder + IMU
        
        for jj = 1:length(rwv_e)
            
            thetaGyroMeas(jj) = Theta0 + randn(1)*(IMUStd); 
            
        end
        
        thetaGyroMeas = thetaGyroMeas';
        
        factor = fEnc/fGPS;
        thetaGyroMeas = decimate(thetaGyroMeas,factor,5);
        
        factor = fEnc/fIMU;
        lwv_e = decimate(lwv_e,factor,5)/dtEnc;
        rwv_e = decimate(rwv_e,factor,5)/dtEnc;
        
        %% EKF

        R = [IMUStd]; 

        V = 1;

        Q = [encStd 0
            0   encStd];

        P_ini = zeros(4,4);
        P = P_ini;
                      
        % Initial Guess of State

        X_hat_prev = [X0, Y0, 0, Theta0]; 

        X_hat_minus_plot = zeros(length(pose_noise),4);

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
            ( -dtEnc/l)     (dtEnc/l)];


        % Measurment Jacobian
        H = [   0 0 0 1   ];


        z = [thetaGyroMeas(i)]; % Use the previous X position
        
        h_0 = H*X_hat_minus; 

       % Kalman Filter Implimentation

       P = A*P*A' + W*Q*W'; 

       K = P*H'*(H*P*H' + V*R*V')^-1;

       %EKF Correction
       X_hat = X_hat_minus  + K*(z - h_0);

       X_kal(i+1,1:4) = X_hat(1:4)';

       P = (eye(4)-K*H)*P; 

       X_hat_prev = X_hat;

    
        end
    
    elseif Option == 3
        
        for jj = 1:length(rwv_e)
            
            thetaGyroMeas(jj) = Theta0 + randn(1)*IMUStd; 
            
        end
        
            thetaGyroMeas = thetaGyroMeas';
        
            factor = fEnc/fGPS;
            thetaGyroMeas = decimate(thetaGyroMeas,factor,5);
        
            factor = fEnc/fIMU;
            lwv_e = decimate(lwv_e,factor,5)/dtEnc;
            rwv_e = decimate(rwv_e,factor,5)/dtEnc;
            
            x_sim_gps = zeros(1,length(XIdeal));
            y_sim_gps = zeros(1,length(XIdeal));

            tt = linspace(0,2*pi,length(XIdeal));

        for q = 1:length(XIdeal)

            x_sim_gps(q) = XIdeal(q) + randn()*(GPSStd) + sinamp*sin(tt(q)/2); 
            y_sim_gps(q) = YIdeal(q) + randn()*(GPSStd);

        end

        factor = fEnc/fGPS;
        x_sim_gps = decimate(x_sim_gps,factor,5);
        y_sim_gps = decimate(y_sim_gps,factor,5);
        
        
        %% EKF

R = [GPSStd 0 0 0 
    0 GPSStd 0 0 
    0 0 GPSStd 0 
    0 0 0 IMUStd]; 
    
V = eye(4);

Q = [encStd 0
    0   encStd];
  
P_ini = zeros(4,4);
P = P_ini;

% Initial Guess of State

X_hat_prev = [X0, Y0, 0, Theta0]; 

X_hat_minus_plot = zeros(length(pose_noise),4);

X_hat_prev(1) = X0;
X_hat_prev(2) = Y0;
X_hat_prev(3) = 0;
X_hat_prev(4) = Theta0;

X_kal = zeros(length(rwv_e),4);

X_kal(1,1:4) = [X_hat_prev(1) X_hat_prev(2) X_hat_prev(3) X_hat_prev(4)];  
    
for i = 1:IMU_points
      
       
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
        ( -dtEnc/l)     (dtEnc/l)];
      
    
    % Measurment Jacobian
    H = [1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1   ];
    
        
    z = [x_sim_gps(i)
        y_sim_gps(i)
        0
        thetaGyroMeas(i)]; % Use the previous X position
         

  
    
    h_0 = H*X_hat_minus; %% from online presentation???
    
   % Kalman Filter Implimentation
     
   P = A*P*A' + W*Q*W'; 
  
   K = P*H'*inv(H*P*H' + V*R*V'); %%%%%
   
   %EKF Correction
   X_hat = X_hat_minus  + K*(z - h_0);
   
   X_kal(i+1,1:4) = X_hat(1:4)';
   
   P = (eye(4)-K*H)*P; 
   
   X_hat_prev = X_hat;
   
end
        
        elseif Option == 4
            
            

            
            
        for jj = 1:length(rwv_e)
            
            thetaGyroMeas(jj) = Theta0 + randn(1)*IMUStd; 
            
        end
        
            thetaGyroMeas = thetaGyroMeas';
            factor = fEnc/fGPS;
            thetaGyroMeas = decimate(thetaGyroMeas,factor,5);
        
            
            
            
        xt = X0;
        yt = Y0;
        theta = Theta0;

        for ii = 1:length(speeds)
            pose_noise(ii,1:4) = [xt, yt, 0, theta];

            xt = pose_noise(ii,1) + lwv_e(i)*dtEnc*cos(pose_noise(ii,4)); 
            yt = pose_noise(ii,2) + lwv_e(i)*dtEnc*sin(pose_noise(ii,4)); 
            theta = pose_noise(ii,4) + (theta_dot_v(rwv_e(ii), lwv_e(ii))); % if both wheels arnt taken into account for dead reckoning, line is straight.

        end
        
        factor = fEnc/fIMU;
        singleEnc = decimate(lwv_e,factor,5)/dtEnc;
        %% EKF

R = [IMUStd]; 
    
V = 1;

Q = [encStd];
  
P_ini = zeros(4,4);
P = P_ini;

% Initial Guess of State

X_hat_prev = [X0, Y0, 0, Theta0]; 

X_hat_minus_plot = zeros(length(pose_noise),4);

X_hat_prev(1) = X0;
X_hat_prev(2) = Y0;
X_hat_prev(3) = 0;
X_hat_prev(4) = Theta0;

X_kal = zeros(length(thetaGyroMeas),4);

X_kal(1,1:4) = [X_hat_prev(1) X_hat_prev(2) X_hat_prev(3) X_hat_prev(4)];  
X_hat_prev_prev = 90*180/pi;    

for i = 2:length(thetaGyroMeas)

    % Assumpition that both wheels do the same thing is too much?
    
    X_hat_minus = [X_hat_prev(1) + singleEnc(i)*dtEnc*cos(X_hat_prev(4))
                   X_hat_prev(2) + singleEnc(i)*dtEnc*sin(X_hat_prev(4))
                   X_hat_prev(3) 
                   X_hat_prev(4) + (thetaGyroMeas(i) - thetaGyroMeas(i-1))  ]; 
       
    
    A = [1  0   0   -(singleEnc(i)*dtEnc)*sin(X_hat_prev(4))
        0   1   0  singleEnc(i)*dtEnc*cos(X_hat_prev(4))
        0   0   1   0
        0   0   0   1]; 
    
    W = [(dtEnc)*cos(X_hat_prev(4))  
        (dtEnc)*sin(X_hat_prev(4)) 
            0 
            0  ]   ;
        
    % Measurment Jacobian
    H = [0 0 0 1];
  
    z = [thetaGyroMeas(i)]; % Use the previous X position

    h_0 = H*X_hat_minus; 
    
   % Kalman Filter Implimentation
     
   P = A*P*A' + W*Q*W'; 
  
   K = P*H'*inv(H*P*H' + V*R*V'); %%%%%
   
   %EKF Correction
   X_hat = X_hat_minus  + K*(z - h_0);
   
   X_kal(i+1,1:4) = X_hat(1:4)';
   
   P = (eye(4)-K*H)*P; 
   
   X_hat_prev = X_hat;
   
end
  
        end
        
if Option == 1
    subplot(2,2,pp)
    hold on

    plot(pose_noise(:,1), pose_noise(:,2),'r')
    scatter(X_kal(:,1),X_kal(:,2),'b','filled')
    plot(X_kal(:,1),X_kal(:,2),'b')
    plot(XIdeal, YIdeal,'k')
    scatter(x_sim_gps,y_sim_gps,'k')
    title(sprintf('2 Encoders + GPS EKF: Sim %s', num2str(pp)),'FontSize', 20)
    xlabel('X Position [mm]')
    ylabel('Y Position [mm]')
    legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','GPS Measurment','Bestoutside', 'Location', 'northeastoutside')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    grid on
    hold off

elseif Option == 2
    subplot(2,2,pp)
    hold on

    plot(pose_noise(:,1), pose_noise(:,2),'r')
    scatter(X_kal(:,1),X_kal(:,2),'b','filled')
    plot(X_kal(:,1),X_kal(:,2),'b')
    plot(XIdeal, YIdeal,'k')
    title(sprintf('2 Encoder + IMU Fusion EKF: Sim %s', num2str(pp)),'FontSize', 20)
    xlabel('X Position [mm]')
    ylabel('Y Position [mm]')
    legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','Bestoutside', 'Location', 'northeastoutside')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    grid on
    hold off

elseif Option == 3
    subplot(2,2,pp)
    hold on

    plot(pose_noise(:,1), pose_noise(:,2),'r')
    scatter(X_kal(:,1),X_kal(:,2),'b','filled')
    plot(X_kal(:,1),X_kal(:,2),'b')
    plot(XIdeal, YIdeal,'k')
    scatter(x_sim_gps,y_sim_gps,'k')
    title(sprintf('2 Encoder + GPS + IMU EKF: Sim %s', num2str(pp)),'FontSize', 20)
    xlabel('X Position [mm]')
    ylabel('Y Position [mm]')
    legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','GPS Measurment','Bestoutside', 'Location', 'northeastoutside')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    grid on
    hold off
    
elseif Option == 4
    subplot(2,2,pp)
    hold on

    plot(pose_noise(:,1), pose_noise(:,2),'r')
    scatter(X_kal(:,1),X_kal(:,2),'b','filled')
    plot(X_kal(:,1),X_kal(:,2),'b')
    plot(XIdeal, YIdeal,'k')
    title(sprintf('1 Encoder + IMU Fusion EKF: Sim %s', num2str(pp)),'FontSize', 20)
    xlabel('X Position [mm]')
    ylabel('Y Position [mm]')
    legend('Dead Reckoning','Kalman Filter Scatter','Kalman Filter Line','Assumed Straight Path','Bestoutside', 'Location', 'northeastoutside')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    grid on
    hold off
    
else 
    
end
    

end






