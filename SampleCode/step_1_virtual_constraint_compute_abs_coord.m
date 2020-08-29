function step_1_virtual_constraint_compute_abs_coord
% This file allows tuning of Bezier coefficients describing the virtual
% constraint of the compass gait walker as per method described in textbook
% by Westervelt (2007).  The Bezier polynomial coefficients generated in
% this file is used in subsequent steps to define the 'nominal' or 'target'
% trajectory.

%% Parameters
g = 9.8;
mh = 10;                %body mass [kg]         worked: 10
m = 5;                  %leg masses [kg]		worked: 5
a = 0.5;                %'shin' length [m]      worked: 0.5
b = 0.5;                %'thigh' length [m]     worked: 0.5
%gamma = 3.1;            %pitch of ramp in deg	worked: 3
gamma = 3;
gamma = deg2rad(gamma); % convert gamma to radians.
l = a + b;              % define length of leg from a and b parameters.

% order of Bez
M = 5;

a_traditional = zeros(1, M+1);

% set up IC / boundary conditions
q1Final = 0.218668812696748; % approx from simulation; end of limit cycle...
q1Init = -0.323388556507581; % because of symmetry and flat ground

period_theta = q1Final - q1Init;
theta_band_between_each_a = period_theta/M;
final_q1_slope = -1.084;
a_traditional(1) = q1Init;
a_traditional(end) = q1Final;

%% HAND TUNE a2/a3/a4 VARIABLES HERE!!
a_traditional(3) = -0.15;    
a_traditional(4) = 0.65;
a_traditional(end-1) = a_traditional(end) - (final_q1_slope*theta_band_between_each_a);  % 0.349679000410283

% theta is the stance leg aka the traditional 'phase' variable
thetaFinal = q1Init; % because of the symmetry of the robot...
thetaInit = q1Final;
q2plus = thetaInit;
q2minus = thetaFinal;

% Various matrices from 
H0 = [1 0];  % what is to be controlled?  controlled_q = H0*q
c = [0 1]; % what is the desired phase variable?  phase_Q = c*q
H = [1 0;0 1];

w = H\[(M/(q2minus-q2plus))*(a_traditional(end) - a_traditional(end-1));1];

% Impact map from textbook
alpha = abs((q1Final - q2minus))/2;  % Calculate alpha, the angle between the legs.

% Angular momentum about point of impact at the instant before impact.
Q_minus = Q_minus_matrix(alpha);
% Angular momentum about point of impact at the instant after impact.
Q_plus = Q_plus_matrix(alpha);

H = (Q_plus)\Q_minus;   % A\B is equiv to inv(A)*B but better

deltaqDot = H;

% Compute a1
a_traditional(2) = (H0*deltaqDot*w*((q2minus-q2plus)/M))/(c*deltaqDot*w) + a_traditional(1);
initial_slope = (a_traditional(2) - a_traditional(1))/theta_band_between_each_a;

a = fliplr(a_traditional) % this is the correct order of bez coef for our phase variable, which is monotonically decreasing

%% plot the guy
figure(); hold on;
func_plot = sprintf('y=-x-2*%f',gamma);
switchSurf_handle = ezplot(func_plot,[-pi/5,pi/5]);

sampleTimes = 100;
phase_all = linspace(thetaInit, thetaFinal, sampleTimes);
h_all = zeros(1,sampleTimes);
for i=1:sampleTimes
    phaseNow = phase_all(i);
    % Remember theta is monotonically decreasing for this model
    phase_var_normalised = (phaseNow-(thetaFinal))/(thetaInit - thetaFinal);
    h_all(i) = bezier(a,phase_var_normalised);
end

virtualCon_handle = plot(phase_all, h_all);


bezCoefPlot = plot(linspace(thetaFinal, thetaInit, M+1),a,'o');

legend([switchSurf_handle,virtualCon_handle,bezCoefPlot],'Switch Surf','Virtual Cons traj','Bezier Poly Coef')
title('Virtual Constraint for actuated nominal trajectory')
save bez_coef

    function sum = bezier(a,s)
        % Compute bezier polynomial output and derivative...
        % Compute the Bezier polynomial
        if isequal(class(a(1)),'sym')
            sum = sym('sum'); sum=sum-sum; % make sure matlab doesn't try to shape things into double
        else
            sum = 0;
        end
        order = length(a) - 1;
        % compute bezier
        for j = 1:(order+1)
            k = j-1;
            sum = sum + a(j)*(s^k)*((1-s)^(order-k))*factorial(order)/(factorial(k)*factorial(order-k));
        end

    end

    function output = Q_minus_matrix(alpha)
        
        output(1,1) = -m*a*b;
        output(1,2) = -m*a*b + (mh*l^2 + 2*m*a*l)*cos(2*alpha);
        output(2,1) = 0;
        output(2,2) = -m*a*b;
    end
%% Q+ matrix
    function output = Q_plus_matrix(alpha)
        
        output(1,1) = m*b*(b - l*cos(2*alpha));
        output(1,2) = m*l*(l - b*cos(2*alpha)) + m*a^2 + mh*l^2;
        output(2,1) = m*b^2;
        output(2,2) = -m*b*l*cos(2*alpha);
    end


end