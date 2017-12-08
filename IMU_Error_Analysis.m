%% Error Analysis
%   How will the position measurment be affected if the gravity vector
%   is measured inacuratly?
%

close all
clear all

r = linspace(0.1,5,10);

g = 9.81;
z_IMU = linspace(9.4,9.81,10); %% z_IMU can not be larger than gravity otherwise... COMPLEX NUMEBRS!
percent_error = 0.01;

errors = zeros(10,2);
theta_deg = zeros(length(z_IMU),1); 
figure(1)
figure(2)
 
for j = 1:length(z_IMU)
    for i = 1:length(r)
        theta = asin(z_IMU(j)/g);
        theta_deg(j) = theta*180/pi;
        
        
        error = g*percent_error;
        theta_wrong = asin((z_IMU(j)-error)/g);

        error =abs(theta-theta_wrong);
        error_deg = error*180/pi;

        circ_error = 2*pi*r(i)*(error_deg/(360));

        errors(i,1) = error;

        errors(i,2) = circ_error;

    end
    
    figure(1)
    hold on
    plot(r, errors(:,2))
    %set(gca,'YTick',0:0.01:1);
    
    grid on
    
    
    title(sprintf('Arc Error with %s%% Error in Z Accel', num2str(percent_error*100)),'FontSize', 18)
    xlabel('Pipe Radius [m]')
    ylabel('Circumference Displacment [m]')
    
    legend(sprintf('\\theta = %s',num2str(theta_deg(1))),sprintf('\\theta = %s',num2str(theta_deg(2))),sprintf('\\theta = %s',num2str(theta_deg(3))),sprintf('\\theta = %s',num2str(theta_deg(4))),sprintf('\\theta = %s',num2str(theta_deg(5))),sprintf('\\theta = %s',num2str(theta_deg(6))),sprintf('\\theta = %s',num2str(theta_deg(7))),sprintf('\\theta = %s',num2str(theta_deg(8))),sprintf('\\theta = %s',num2str(theta_deg(9))),sprintf('\\theta = %s',num2str(theta_deg(10))),'Location', 'northeastoutside')   
    hold off
    
   
    
end

r = linspace(0.1,1,10);

g = 9.81;
z_IMU = linspace(0,9.81,2000); 
percent_error = 0.001;

errors = zeros(10,2);
theta = zeros(length(z_IMU),1);
theta_deg = zeros(length(z_IMU),1); 
percentage_error = zeros(length(z_IMU),1);    

for j = 1:length(r)
    
    for i = 1:length(z_IMU)
        theta(i) = asin(z_IMU(i)/g);
        theta_deg(i) = theta(i)*180/pi;
        
        
        error = g*percent_error;
        theta_wrong = asin((z_IMU(i)-error)/g);

        error =abs(theta(i)-theta_wrong);
        error_deg = error*180/pi;

        circ_error = 2*pi*r(j)*(error_deg/(360));

        errors(i,1) = error;

        errors(i,2) = circ_error;
        
        percentage_error(i) = circ_error/(2*pi*r(j));
        
    end
    
    figure(2)
    hold on
    plot(theta_deg, errors(:,2))
    set(gca,'YTick',0:0.005:1);
    
    grid on
        
    title(sprintf('Arc Error with %s%% Error in Z Accel', num2str(percent_error*100)),'FontSize', 18)
    xlabel('\theta [°]')
    ylabel('Circumference Displacment [m]')
    xlim([0 90])
    
    legend(sprintf('r = %s m',num2str(r(1))),sprintf('r = %s m',num2str(r(2))),sprintf('r = %s m',num2str(r(3))),sprintf('r = %s m',num2str(r(4))),sprintf('r = %s m',num2str(r(5))),sprintf('r = %s m',num2str(r(6))),sprintf('r = %s m',num2str(r(7))),sprintf('r = %s m',num2str(r(8))),sprintf('r = %s m',num2str(r(9))),sprintf('r = %s m',num2str(r(10))),'Location', 'northwest')   
    hold off
    
    figure(3)
    
    hold on
    plot(theta_deg, percentage_error)
    %set(gca,'YTick',0:0.01:1);
    
    grid on

    title(sprintf('Error as %% of Circumference\n with %s%% Error in Z Accel', num2str(percent_error*100)),'FontSize', 18)
    xlabel('\theta [°]')
    ylabel('%% Error')
    xlim([0 90])
    
    legend(sprintf('r = %s m',num2str(r(1))),sprintf('r = %s m',num2str(r(2))),sprintf('r = %s m',num2str(r(3))),sprintf('r = %s m',num2str(r(4))),sprintf('r = %s m',num2str(r(5))),sprintf('r = %s m',num2str(r(6))),sprintf('r = %s m',num2str(r(7))),sprintf('r = %s m',num2str(r(8))),sprintf('r = %s m',num2str(r(9))),sprintf('r = %s m',num2str(r(10))),'Location', 'northwest')   
    hold off
    
    
end













