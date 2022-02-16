%% Problem 1
METERS_PER_KM = 1000;
SECONDS_PER_HOUR = 60*60;
Cd = 0.32; % Aerodynamic drag coefficient
A = 1.71; % Area
air_density = 1.293;
g = 9.8; % gravity
m0 = 1100; % kg
m1 = 150; % kg
up = 0.013; % rolling coefficient
grade = 0.05; % Grade percentage
alpha = atan(grade); % Radians
mass_total = m0 + m1;
velocity = 100*METERS_PER_KM/SECONDS_PER_HOUR;
% Weight
W = mass_total*g;
% Gradient resistance
Fg = W*sin(alpha);
% Aerodynamic drag force
Fa = Cd*A*air_density*(velocity^2)/2;
try
    % Tire rolling resistance
    Ff = up*W*(1-(up*W/(4*Ct*St)));
catch
    % Simplified tire rolling resistance formula if parameters are missing
    % from above equation
    Ff = up*W;
end
Ft = Fa + Ff + Fg;
fprintf('Problem 1:\n');
fprintf('Aerodynamic drag force (Fa) = %.2f N\n', Fa);
fprintf('Tire rolling resistance (Ff) = %.2f N\n', Ff);
fprintf('Gradient resistance (Fg) = %.2f N\n', Fg);
fprintf('Total resistance force: Fa+Ff+Fg = %.02f KN\n', Ft/1000);
fprintf('--------------------------------------------------\n');

%% Problem 2
close all;
xMax = 5;
yMax = 3.4;
theta = 23;
fprintf('Problem 2:\n');
fprintf('Friction ellipse equation: (Fx/Fxmax)^2 + (Fy/Fymax)^2 = 1\n');
displayEllipse(xMax, yMax, theta)
xMax = 12; yMax = 6; theta = 45;
displayEllipse(xMax, yMax, theta);
fprintf('--------------------------------------------------\n');

%% Problem 3
L = 4; % Meters
R = 400; % Meters
alpha_f = deg2rad(5); % Converting to radians
alpha_r = deg2rad(3); % Converting to radians

Steering_Angle_High_Speed = L/R + alpha_f - alpha_r;
Steering_Angle_Low_Speed = L/R; % Assuming minimal tire slip angles
fprintf('Problem 3:\n');
fprintf('1) Steering Angle at high speed: %.1f degrees\n',rad2deg(Steering_Angle_High_Speed));
fprintf('2) Steering Angle at low speed: %.1f degrees\n',rad2deg(Steering_Angle_Low_Speed));
fprintf('--------------------------------------------------\n');

%% Problem 4
L = 4;
Wf = 2500; % Pounds
Wr = 1500; % Pounds
Caf = 200; % Pounds
Car = 150; % Pounds

Kus = Wf/Caf - Wr/Car;
if Kus > 0
    Steering_Condition = 'Understeer';
elseif Kus < 0
    Steering_Condition = 'Oversteer';
else
    Steering_Condition = 'Neutral steer';
end

Vcrit = sqrt(g*L/abs(Kus));
fprintf('Problem 4:\n');
fprintf('1) Understeer coefficient: %.02f\n',Kus); 
fprintf('2) The vehicle state is %s\n', Steering_Condition);
fprintf('3) The critical speed is: %.02f m/s\n', Vcrit);

%% Problem 2 Ellipse function
function [x, y] = getEllipse(a, b)
    % Plot range
    N = 100;
    theta = 0:1/N:2*pi+1/N;

    % Equation of ellipse in Parametric Form 
    % https://testbook.com/learn/maths-equation-of-ellipse/#:~:text=Equation%20of%20Ellipse%20Solved%20Examples%201%20Solution%3A%20The,both%20the%20foci%20rest%20on%20the%20Y-axis.%20
    x = a*cos(theta); 
    y = b*sin(theta);
end

%% Problem 2 Ellipse plot
function[] = displayEllipse(xMax, yMax, theta)
    figure();
    h = zeros(1,2);
    thetaRad = deg2rad(theta);
    [x,y] = getEllipse(xMax, yMax);
    h(1) = plot(x,y);
    hold on
    
    % Fy, Fx line
    xPoint = xMax*cos(thetaRad); 
    yPoint = yMax*sin(thetaRad);
    h(2) = plot([0 xPoint], [0 yPoint]);
    hold on
    
    % Create grid
    plot([-xMax xMax], [0 0],'k');
    plot([0 0], [-yMax yMax], 'k');
    hold on
    plot([0 xPoint], [yPoint yPoint],'k', 'LineStyle', '--');
    plot([xPoint xPoint], [0 yPoint],'k', 'LineStyle', '--');
    txt = 'Fy max';
    text(0,(yMax+yMax/10),txt,'HorizontalAlignment','center');
    txt = 'Fx max';
    text(xMax,0,txt,'HorizontalAlignment','left');
    txt = 'Fy \rightarrow';
    text(0,yPoint,txt,'HorizontalAlignment','right');
    txt = 'Fx \uparrow     ';
    text(xPoint,-0.25,txt,'HorizontalAlignment','center');
    
    % Make the x,y axis the same size
    m = max(xMax, yMax);
    axis([-m m -m m]);
    title('Friction ellipse');
    xlabel('Braking or tractive force (KN)');
    ylabel('Cornering force (KN)')
    grid;
    legend(h(1:2),'Friction ellipse', 'Fx, Fy vector');
    
    % Return Fx, Fy
    fprintf('Given FxMax = %.1f KN, FyMax = %.1f KN, and theta = %i degrees, Fx = %.2f KN and Fy = %.2f KN\n',xMax, yMax, theta, xPoint, yPoint);
end