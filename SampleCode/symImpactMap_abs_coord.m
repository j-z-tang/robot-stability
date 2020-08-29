function [deltaq,deltaqDot] = symImpactMap_abs_coord(q,p)

if nargin<2
    % default values
    % Robot parameters as per Westervelt textbk p.65
    g = 9.8;
    mh = 10;                %body mass [kg]         worked: 10
    m = 5;                  %leg masses [kg]		worked: 5
    a = 0.5;                %'shin' length [m]      worked: 0.5
    b = 0.5;                %'thigh' length [m]     worked: 0.5
    %gamma = 3.1;            %pitch of ramp in deg	worked: 3
    gamma = 3;
    gamma = deg2rad(gamma); % convert gamma to radians.
    l = a + b;              % define length of leg from a and b parameters.
    
else
    l= p.l;
    a = p.a_len;
    b = p.b_len ;
    mh = p.mh ;
    m = p.m_leg;
    g = p.g0;
    gamma = p.gamma ;
end
% Calculate alpha, the angle between the legs.
alpha = ((q(1) - q(2)))/2;

% if ~isequal(class(q),'sym')
%     if alpha>0
%         warning('alpha>0')
%     end
% end


% Computes impact map for symmetrical robot
deltaq = [0 1; 1 0];

% Using conservation of angular momentum, the state before and after
% impact can be calculated.

% Angular momentum about point of impact at the instant before impact.
Q_minus = Q_minus_matrix(alpha);
% Angular momentum about point of impact at the instant after impact.
Q_plus = Q_plus_matrix(alpha);

deltaqDot = (Q_plus)\Q_minus;   % A\B is equiv to inv(A)*B but better

%% Q- matrix
    function output = Q_minus_matrix(alpha)
        if isequal(class(alpha),'sym')
            output = sym('Q',[2,2]);
        end
        output(1,1) = -m*a*b;
        output(1,2) = -m*a*b + (mh*l^2 + 2*m*a*l)*cos(2*alpha);
        output(2,1) = 0;
        output(2,2) = -m*a*b;
    end
%% Q+ matrix
    function output = Q_plus_matrix(alpha)
        if isequal(class(alpha),'sym')
            output = sym('Q',[2,2]);
        end
        output(1,1) = m*b*(b - l*cos(2*alpha));
        output(1,2) = m*l*(l - b*cos(2*alpha)) + m*a^2 + mh*l^2;
        output(2,1) = m*b^2;
        output(2,2) = -m*b*l*cos(2*alpha);
    end

end