function step_0_compute_model_abs_coord
% This file generates the required functions encoding the dynamics of the
% compass gait walker.

x = sym('x',[4,1]);
x = sym(x,'real');
u = sym('u','real');
length = sym('length','real');
a_len = sym('a_len','real');
b_len =sym('b_len','real');
mh = sym('mh','real');
m_leg = sym('m_leg','real');
g0 =sym('g0','real');
gamma = sym('gamma','real');
phaseInit = sym('phaseInit');
phaseFinal = sym('phaseFinal');
a0= sym('a0');
a1= sym('a1');
a2= sym('a2');
a3= sym('a3');
a4= sym('a4');
a5= sym('a5');
a = [a0,a1,a2,a3,a4,a5];
P_gain = sym('P_gain');
D_gain = sym('D_gain');

%% control
display('Computing control input for full f')
input_PD_FL_abs_full_f = simplify(control_PD_FL_absCoord_full_f(x));

%% compute model with u
display('Computing dynamics for full f')
f_full_abs=  simplify(EoM_full_absCoord([],x, u));

dfdx_full_abs = simplify(jacobian(f_full_abs, x));
dfdu_full_abs = simplify( jacobian(f_full_abs, u) );

%% write to file...
param_list = {'g0','p.g0';
    'length','p.l';
    'a_len','p.a_len';
    'b_len','p.b_len';
    'mh','p.mh';
    'm_leg','p.m_leg';
    'gamma','p.gamma';
    'phaseInit','p.phaseInit';
    'phaseFinal','p.phaseFinal';
    'a0','p.A_vec(1)';
    'a1','p.A_vec(2)';
    'a2','p.A_vec(3)';
    'a3','p.A_vec(4)';
    'a4','p.A_vec(5)';
    'a5','p.A_vec(6)';
    'P_gain','p.P_gain';
    'D_gain','p.D_gain'};

list_x = {'x1','x(1)'; 'x2','x(2)'; 'x3','x(3)'; 'x4','x(4)'};

list_u = {'u','u'};

write_fcn('f_full_abs.m',{'x','u','p'},[list_x; list_u; param_list],{f_full_abs,'f_full_abs'});

write_fcn('dfdx_full_abs.m',{'x','u','p'},[list_x; list_u; param_list],{dfdx_full_abs,'dfdx_full_abs'});
write_fcn('dfdu_full_abs.m',{'x','p'},[list_x; param_list],{dfdu_full_abs,'dfdu_full_abs'});

write_fcn('input_PD_FL_abs_full_f.m',{'x','p'},[list_x; param_list],{input_PD_FL_abs_full_f,'input_PD_FL_abs_full_f'});

%% Equation of Motion for compasss....
    function output = EoM_full_absCoord(~,q, u)
        
        % Rename state variables for clarity.
        theta(1,1) = q(1);  % theta swing
        theta(2,1) = q(2);  % theta stance
        theta_dot(1,1) = q(3);
        theta_dot(2,1) = q(4);
        
        % Calculate the mass, coriolis and gravity matrices.
        H = D_Matrix(theta);
        C = C_Matrix(theta,theta_dot);
        G = G_Matrix(theta);
        
        B = [1; -1];
        
        % Calculate the derivative of the state variable.
        q_dot(1:2) = theta_dot;
        
        q_dot(3:4) =  -H\(C*theta_dot + G - B*u); % with input signal
        
        % Output of the function is the derivative.
        output = q_dot';
    end

    function [control, y] = control_PD_FL_absCoord_full_f(q)
        
        q1=q(1); q2=q(2); q3=q(3); q4=q(4);
        
        phaseVar = q(2);
        % theta is sweeping the other way!
        phase_var_normalised = (phaseVar-(phaseFinal))/(phaseInit - phaseFinal);
        y = q(1) - bezier(a,phase_var_normalised);
        
        % output function (virtual constraint)
        Lfh = q3 + q4*((5*a0*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^4)/(phaseInit - phaseFinal) - (5*a1*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^4)/(phaseInit - phaseFinal) + (5*a4*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (5*a5*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (20*a1*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^3*(phaseFinal - q2))/(phaseInit - phaseFinal)^2 - (20*a3*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^4 + (20*a4*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^4 + (10*a2*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^3*(2*phaseFinal - 2*q2))/(phaseInit - phaseFinal)^2 + (30*a2*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^3 - (30*a3*((phaseFinal - q2)/(phaseInit - phaseFinal) + 1)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^3);
        
        LgLfh = ((a_len*cos(q1 - q2) - b_len + b_len*cos(q1 - q2))*((5*a0*(phaseInit - q2)^4)/(phaseInit - phaseFinal)^5 - (5*a1*(phaseInit - q2)^4)/(phaseInit - phaseFinal)^5 + (5*a4*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (5*a5*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (20*a1*(phaseInit - q2)^3*(phaseFinal - q2))/(phaseInit - phaseFinal)^5 - (20*a3*(phaseInit - q2)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^5 + (20*a4*(phaseInit - q2)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^5 + (10*a2*(phaseInit - q2)^3*(2*phaseFinal - 2*q2))/(phaseInit - phaseFinal)^5 + (30*a2*(phaseInit - q2)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^5 - (30*a3*(phaseInit - q2)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^5))/(b_len*(2*a_len^2*m_leg + a_len^2*mh + b_len^2*m_leg + b_len^2*mh + 2*a_len*b_len*m_leg + 2*a_len*b_len*mh - a_len^2*m_leg*cos(q1 - q2)^2 - b_len^2*m_leg*cos(q1 - q2)^2 - 2*a_len*b_len*m_leg*cos(q1 - q2)^2)) + (2*a_len^2*m_leg + a_len^2*mh + b_len^2*m_leg + b_len^2*mh + 2*a_len*b_len*m_leg + 2*a_len*b_len*mh - b_len^2*m_leg*cos(q1 - q2) - a_len*b_len*m_leg*cos(q1 - q2))/(b_len^2*m_leg*(2*a_len^2*m_leg + a_len^2*mh + b_len^2*m_leg + b_len^2*mh + 2*a_len*b_len*m_leg + 2*a_len*b_len*mh - a_len^2*m_leg*cos(q1 - q2)^2 - b_len^2*m_leg*cos(q1 - q2)^2 - 2*a_len*b_len*m_leg*cos(q1 - q2)^2));
        
        LfLfh = (((5*a0*(phaseInit - q2)^4)/(phaseInit - phaseFinal)^5 - (5*a1*(phaseInit - q2)^4)/(phaseInit - phaseFinal)^5 + (5*a4*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (5*a5*(phaseFinal - q2)^4)/(phaseInit - phaseFinal)^5 - (20*a1*(phaseInit - q2)^3*(phaseFinal - q2))/(phaseInit - phaseFinal)^5 - (20*a3*(phaseInit - q2)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^5 + (20*a4*(phaseInit - q2)*(phaseFinal - q2)^3)/(phaseInit - phaseFinal)^5 + (10*a2*(phaseInit - q2)^3*(2*phaseFinal - 2*q2))/(phaseInit - phaseFinal)^5 + (30*a2*(phaseInit - q2)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^5 - (30*a3*(phaseInit - q2)^2*(phaseFinal - q2)^2)/(phaseInit - phaseFinal)^5)*(3*a_len*g0*m_leg*sin(q2) - b_len*g0*m_leg*sin(2*q1 - q2) - 2*b_len^2*m_leg*q3^2*sin(q1 - q2) - a_len*g0*m_leg*sin(2*q1 - q2) + 2*a_len*g0*mh*sin(q2) + b_len*g0*m_leg*sin(q2) + 2*b_len*g0*mh*sin(q2) + a_len^2*m_leg*q4^2*sin(2*q1 - 2*q2) + b_len^2*m_leg*q4^2*sin(2*q1 - 2*q2) - 2*a_len*b_len*m_leg*q3^2*sin(q1 - q2) + 2*a_len*b_len*m_leg*q4^2*sin(2*q1 - 2*q2)))/(3*a_len^2*m_leg + 2*a_len^2*mh + b_len^2*m_leg + 2*b_len^2*mh + 2*a_len*b_len*m_leg + 4*a_len*b_len*mh - a_len^2*m_leg*cos(2*q1 - 2*q2) - b_len^2*m_leg*cos(2*q1 - 2*q2) - 2*a_len*b_len*m_leg*cos(2*q1 - 2*q2)) - (20*q4^2*(a0*phaseInit^3 - 2*a1*phaseInit^3 + a2*phaseInit^3 - a3*phaseFinal^3 + 2*a4*phaseFinal^3 - a5*phaseFinal^3 - a0*q2^3 + 5*a1*q2^3 - 10*a2*q2^3 + 10*a3*q2^3 - 5*a4*q2^3 + a5*q2^3 - 3*a1*phaseInit^2*phaseFinal + 3*a2*phaseInit*phaseFinal^2 + 6*a2*phaseInit^2*phaseFinal - 6*a3*phaseInit*phaseFinal^2 - 3*a3*phaseInit^2*phaseFinal + 3*a4*phaseInit*phaseFinal^2 + 3*a0*phaseInit*q2^2 - 3*a0*phaseInit^2*q2 - 12*a1*phaseInit*q2^2 + 9*a1*phaseInit^2*q2 + 18*a2*phaseInit*q2^2 - 9*a2*phaseInit^2*q2 - 12*a3*phaseInit*q2^2 + 3*a3*phaseInit^2*q2 + 3*a4*phaseInit*q2^2 - 3*a1*phaseFinal*q2^2 + 12*a2*phaseFinal*q2^2 - 3*a2*phaseFinal^2*q2 - 18*a3*phaseFinal*q2^2 + 9*a3*phaseFinal^2*q2 + 12*a4*phaseFinal*q2^2 - 9*a4*phaseFinal^2*q2 - 3*a5*phaseFinal*q2^2 + 3*a5*phaseFinal^2*q2 + 6*a1*phaseInit*phaseFinal*q2 - 18*a2*phaseInit*phaseFinal*q2 + 18*a3*phaseInit*phaseFinal*q2 - 6*a4*phaseInit*phaseFinal*q2))/(phaseInit - phaseFinal)^5 - (2*a_len^2*g0*m_leg*sin(q1 - 2*q2) + a_len^2*g0*mh*sin(q1 - 2*q2) + b_len^2*g0*m_leg*sin(q1 - 2*q2) + b_len^2*g0*mh*sin(q1 - 2*q2) - 4*a_len^3*m_leg*q4^2*sin(q1 - q2) - 2*a_len^3*mh*q4^2*sin(q1 - q2) - 2*b_len^3*m_leg*q4^2*sin(q1 - q2) - 2*b_len^3*mh*q4^2*sin(q1 - q2) + b_len^3*m_leg*q3^2*sin(2*q1 - 2*q2) + 2*a_len^2*g0*m_leg*sin(q1) + a_len^2*g0*mh*sin(q1) + b_len^2*g0*m_leg*sin(q1) + b_len^2*g0*mh*sin(q1) + 3*a_len*b_len*g0*m_leg*sin(q1 - 2*q2) + 2*a_len*b_len*g0*mh*sin(q1 - 2*q2) + 2*a_len*b_len^2*m_leg*q3^2*sin(2*q1 - 2*q2) + a_len^2*b_len*m_leg*q3^2*sin(2*q1 - 2*q2) - 6*a_len*b_len^2*m_leg*q4^2*sin(q1 - q2) - 8*a_len^2*b_len*m_leg*q4^2*sin(q1 - q2) - 6*a_len*b_len^2*mh*q4^2*sin(q1 - q2) - 6*a_len^2*b_len*mh*q4^2*sin(q1 - q2) + a_len*b_len*g0*m_leg*sin(q1) + 2*a_len*b_len*g0*mh*sin(q1))/(b_len*(3*a_len^2*m_leg + 2*a_len^2*mh + b_len^2*m_leg + 2*b_len^2*mh + 2*a_len*b_len*m_leg + 4*a_len*b_len*mh - a_len^2*m_leg*cos(2*q1 - 2*q2) - b_len^2*m_leg*cos(2*q1 - 2*q2) - 2*a_len*b_len*m_leg*cos(2*q1 - 2*q2)));
        
        PD_term = -P_gain*y - D_gain*Lfh;
        
        control = LgLfh\(PD_term-LfLfh);
        
    end

    function sum = bezier(a,s)
        % Compute bezier polynomial output and derivative...
        if isequal(class(a(1)),'sym')
            sum = sym('sum'); sum=sum-sum; % make sure matlab doesn't try to shape things into double
        else
            sum = 0;
        end
        order = max(size(a)) - 1;
        % compute bezier
        for j = 1:(order+1)
            k = j-1;
            sum = sum + a(j)*(s^k)*((1-s)^(order-k))*factorial(order)/(factorial(k)*factorial(order-k));
        end
    end


    function v = PFL_feedback_control(q,nominal_u)
        theta(1,1) = q(1);  % theta stance
        theta(2,1) = q(2);  % theta swing
        theta_dot(1,1) = q(3); % stance
        theta_dot(2,1) = q(4); % theta swing
        
        H = D_Matrix(theta);
        C = C_Matrix(theta,theta_dot);
        G = G_Matrix(theta);
        
        h_vec = C*theta_dot;
        
        H_bar_11 = H(1,1) - H(1,2)*(H(2,2)\H(2,1));
        h_vec_bar_1 = h_vec(1) - H(1,2)*(H(2,2)\h_vec(2));
        phi_bar_1 = G(1) - H(1,2)*(H(2,2)\G(2));
        
        v = H_bar_11 \ (  (1+(H(1,2)/H(2,2)) )*nominal_u - h_vec_bar_1 - phi_bar_1);
    end

    function output = EoM_PFL_absCoord(~,q,lin_v)
        
        % Rename state variables for clarity.
        theta(1,1) = q(1);  % theta stance
        theta(2,1) = q(2);  % theta swing
        theta_dot(1,1) = q(3); % stance
        theta_dot(2,1) = q(4); % theta swing
        %         theta_dot(3,1) = q(6);
        
        % Calculate the mass, coriolis and gravity matrices.
        H = D_Matrix(theta);
        C = C_Matrix(theta,theta_dot);
        G = G_Matrix(theta);
        
        h_vec = C*theta_dot;
        
        H_bar_11 = H(1,1) - H(1,2)*(H(2,2)\H(2,1));
        h_vec_bar_1 = h_vec(1) - H(1,2)*(H(2,2)\h_vec(2));
        phi_bar_1 = G(1) - H(1,2)*(H(2,2)\G(2));
        
        q_dot(1:2) = theta_dot;
        q_dot(3) = lin_v; %q1DotDot
        q_dot(4) = -H(2,2) \ ( H(2,1)*lin_v + h_vec(2) + G(2)  +   (1+H(1,2)/H(2,2))\(H_bar_11*lin_v + h_vec_bar_1 + phi_bar_1) );
        
        % Output of the function is the derivative.
        output = q_dot';
    end



%% Dynamics help function
    function output=D_Matrix(theta)
        if isequal(class(theta),'sym')
            output = sym('D',[2,2]);
        end
        output(1,1) = m_leg*b_len^2;
        output(1,2) = -m_leg*length*b_len*cos(theta(2)-theta(1));
        output(2,1) = -m_leg*length*b_len*cos(theta(2)-theta(1));
        output(2,2) = (mh + m_leg)*length^2 + m_leg*a_len^2;
    end

    function output=C_Matrix(theta,theta_dot)
        if isequal(class(theta),'sym')
            output = sym('D',[2,2]);
        end
        
        output(1,1) = 0;
        output(1,2) = m_leg*length*b_len*sin(theta(2)-theta(1))*theta_dot(2);
        output(2,1) = -m_leg*length*b_len*sin(theta(2)-theta(1))*theta_dot(1);
        output(2,2) = 0;
    end

    function output = G_Matrix(theta)
        % theta(1)=swing angle
        % theta(2)=stance angle
        if isequal(class(theta),'sym')
            output = sym('D',[2,1]);
        end
        output(1,1) = m_leg*b_len*g0*sin(theta(1));
        output(2,1) = -(mh*length + m_leg*a_len + m_leg*length)*g0*sin(theta(2));
    end


end