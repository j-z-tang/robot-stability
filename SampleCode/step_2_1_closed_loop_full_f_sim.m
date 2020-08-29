function step_2_1_closed_loop_full_f_sim
% simulate closed loop behaviour of compass walker using virtual constraints feedback
% controller outlined in Westervelt textbook.  This generates the 'open
% loop' trajectory and control input around which the transverse
% coordinates will be constructed.

animate = false;

sim_steps = 60; % number of steps simulated

saveFolder = 'alpha_2_1_closed_loop_full_abs_coord_unactuated';

if ~exist(saveFolder,'dir')
    mkdir(saveFolder)
end

coordTypeAll = {'abs','rel'};
coordType = coordTypeAll{1};

% the user can select the following to decide whether the nominal
% 'target' trajectory is actuated (using virtual constraint designed in
% Step 1 of the code, or use the natural unactuated trajectory as the
% nominal 'target trajectory.
controllerTypeAll = {'f_full_feedback','no_input'};
controllerType = controllerTypeAll{2};

update_file = true;

loadedVar = load('bez_coef.mat','a','thetaInit','thetaFinal'); % bezier poly describing virtual constraint
A_vec = loadedVar.a;
thetaInit = loadedVar.thetaInit;
thetaFinal = loadedVar.thetaFinal;

%% control parameters to find nominal 'open loop' input for given trajectory
% for finite-time control
control_ep = 0.1;
control_alpha = 0.9;
P_gain = 500; 
D_gain = 20;

% desired q3 angle...
q3d = pi/4;
phaseInit = thetaInit;
phaseFinal = thetaFinal;

% Robot parameters as per Westervelt textbk p.65
switch coordType
    case 'rel'
        l = 1;
        l_c = 0.8;
        I_c = 0.03;
        m = 0.3;
        g0 = 9.81;
        slope_alpha = deg2rad(3);
        
        p.l = l;
        p.l_c = l_c;
        p.I_c = I_c;
        p.m = m;
        p.g0 = g0;
        p.slope_alpha = slope_alpha;
    case 'abs'
        
        g = 9.8;
        mh = 10;                %body mass [kg]         worked: 10
        m = 5;                  %leg masses [kg]		worked: 5
        a = 0.5;                %'shin' length [m]      worked: 0.5
        b = 0.5;                %'thigh' length [m]     worked: 0.5
        %gamma = 3.1;            %pitch of ramp in deg	worked: 3
        gamma = 3;
        gamma = deg2rad(gamma); % convert gamma to radians.
        l = a + b;              % define length of leg from a and b parameters.
        
        p.l = l;
        p.a_len = a;
        p.b_len = b;
        p.mh = mh;
        p.m_leg = m;
        p.g0 =g;
        p.gamma = gamma;
end

p.A_vec = A_vec;
p.phaseInit = phaseInit;
p.phaseFinal = thetaFinal;
p.P_gain = P_gain;
p.D_gain = D_gain;

switch coordType
    case 'abs'; 
        
        symImpactMap = @(x,p) symImpactMap_abs_coord(x,p);
        foot_impact = @(t,x) foot_impact_abs([],x,p);
    case 'rel'; 
        
        symImpactMap = @(x,p) symImpactMap_rel_coord(x,p);
        foot_impact = @(t,x) foot_impact_rel([],x);
end


%% select the right dynamics
% inline funcs
switch controllerType
    case 'f_full_feedback'
        switch coordType
            case 'rel'
                feedbackLaw = @(x,param) input_PD_FL_relCoord_full_f(x,param);
                f_full = @(x,u) f_full_rel(x,u,p);
                dfdx_full = @(x,u) dfdx_full_rel(x,u,p);
                dfdu_full = @(x) dfdu_full_rel(x,p);
            case 'abs'
                feedbackLaw = @(x,param) input_PD_FL_abs_full_f(x,param);
                f_full = @(x,u) f_full_abs(x,u,p);
                dfdx_full = @(x,u) dfdx_full_abs(x,u,p);
                dfdu_full = @(x) dfdu_full_abs(x,p);
        end
        
    case 'no_input'
        switch coordType
            case 'rel'
                feedbackLaw = @(x,param) zeros(1,1);
                f_full = @(x,u) f_full_rel(x,u,p);
                dfdx_full = @(x,u) dfdx_full_rel(x,u,p);
                dfdu_full = @(x) dfdu_full_rel(x,p);
            case 'abs'
                feedbackLaw = @(x,param) zeros(1,1);
                f_full = @(x,u) f_full_abs(x,u,p);
                dfdx_full = @(x,u) dfdx_full_abs(x,u,p);
                dfdu_full = @(x) dfdu_full_abs(x,p);
        end
end

%% simulate full dynamics
phase1 = figure(); hold on; 
plot_dims = [2,1]; % which three dimensions of state space to plot


ode45_max_time = 1;
switch coordType
    case 'rel'
        initial_conditions = [phaseInit*2, phaseInit, 1, 1.09];  % WORKS!!! for rel_coord_with_input
    case 'abs'
        initial_conditions = [ phaseFinal    phaseInit   -0.3772 ,  -1.0918]; % works for abs coord
end

% plot initial condition point
plot(initial_conditions(:,plot_dims(1)), initial_conditions(:,plot_dims(2)),'o');
title('Phase portrait of trajectory')
xlabel('$\theta_1$','Interpreter','LaTex'); ylabel('$\theta_2$','Interpreter','LaTex');

odeParam = odeset('Events',foot_impact, 'RelTol',1e-9); % termination condition/settings

CompassStep = struct('t',cell(1,sim_steps),'q',cell(1,sim_steps));
control_input_record = cell(sim_steps);

for run=1:sim_steps
    
    % swing dynamics
    [CompassStep(run).t, CompassStep(run).q] = ode45(@(t,x) f_full(x,feedbackLaw(x,p)), [0,ode45_max_time], initial_conditions, odeParam);
    figure(phase1); plot(CompassStep(run).q(:,plot_dims(1)), CompassStep(run).q(:,plot_dims(2)),'b');

    % switching (impact)
    q_final = CompassStep(run).q(end,:);
    
    [deltaq, deltaqDot] = symImpactMap(q_final,p);
    
    initial_conditions = [deltaq*[q_final(1:2)']; deltaqDot*[q_final(3:4)]']';
    
    switchLine = [q_final;initial_conditions];
    figure(phase1); plot(switchLine(:,plot_dims(1)), switchLine(:,plot_dims(2)),'b--');
end


%% re-run limit cycle to so we have equally spaced time...
limCycleLength = CompassStep(end).t(end);
sampleTimes = 2000; % number of samples points within a step

limCycle = struct();

[limCycle.t, limCycle.q] = ode45(@(t,x) f_full(x,feedbackLaw(x,p)), linspace(0,limCycleLength,sampleTimes), initial_conditions, odeParam);


%% nominal control action and store u_star
control_record=zeros(1,length(CompassStep(run).t));
y_record=zeros(1,length(CompassStep(run).t));

for i= 1:length(limCycle.t)
%      [control_record(i), y_record(i)] = control_simple(CompassStep(run).q(i,:));
    control_record(i) = feedbackLaw(limCycle.q(i,:),p);
end
nominal_u_spline = spline(limCycle.t, control_record);
nominal_u_pp = @(t) ppval(nominal_u_spline,t);

y_record;

%% store x_star
xStarSpline = spline(limCycle.t, limCycle.q');
xDesiredpp = @(t) ppval(xStarSpline,t);

switch coordType
    case 'rel'
        xDotDesiredpp = @(t) f_full(xDesiredpp(t),input_PD_FL_relCoord_full_f(xDesiredpp(t),p));
    case 'abs'
        xDotDesiredpp = @(t) f_full(xDesiredpp(t),input_PD_FL_abs_full_f(xDesiredpp(t),p));
end


% plot in absolute coordinates...
if coordType ~= 'abs'
swingAbsCoord = limCycle.q(:,2)-limCycle.q(:,1); % stance - swing
figure(); hold on;
ezplot('y=-x',[-pi/5,pi/5]);
plot(limCycle.q(:,2), swingAbsCoord,'r');
title('Abs Coordinate');
end

if update_file
    clear phase1
    save(sprintf('%s/result_0_2_closed_loop_sim_%s.mat',saveFolder,controllerType));
end

%% animation (activate at the top)
if animate
    
    l = 1;
    L = 0.5;
    
    animationWrtTime = figure('units','normalized','PaperType','a3');
    hold on
    axis equal
    axis([-l-.5 l+.5 -.5 L+l+.5]);
    
    hHip = plot(10,10,'ko','MarkerSize',10,'MarkerFaceColor','g');
    hStance = plot([10 11],[10 11],'r-','LineWidth',2);
    hSwing = plot([10 11],[10 11],'b-','LineWidth',2);
%     hTorso = plot([10 11],[10 11],'k-','LineWidth',2);
    hGround = plot([-l-.5 l+.5],[0 0],'k');
    
    for run=1:sim_steps
        for timeIndex = 1:length(CompassStep(run).t)
            
            timeNow = CompassStep(run).t(timeIndex);
            
            % Just begin with the hip position
            qStanceNow = CompassStep(run).q(timeIndex,2);
            qSwingNow = CompassStep(run).q(timeIndex,1);
            
            hipPosGround = [0; cos(qStanceNow)*(l)];
            % keep stance leg position at y axis
            stanceLegEndPos = [-sin(qStanceNow)*(l); 0];
            swingLegEndPos = [hipPosGround(1)-sin(-qSwingNow+qStanceNow)*(l); hipPosGround(2)-cos(-qSwingNow+qStanceNow)*(l)];

            % update robot visualization
            figure(animationWrtTime)
            set(hStance,'xdata',[stanceLegEndPos(1) hipPosGround(1)],'ydata',[stanceLegEndPos(2) hipPosGround(2)])
            set(hSwing,'xdata',[swingLegEndPos(1) hipPosGround(1)],'ydata',[swingLegEndPos(2) hipPosGround(2)])
%             set(hTorso,'xdata',[torsoEndPos(1) hipPosGround(1)],'ydata',[torsoEndPos(2) hipPosGround(2)])
            set(hHip,'xdata',hipPosGround(1),'ydata',hipPosGround(2))
            title(sprintf('t = %f s',timeNow));
            %     title(['Step:' num2str(k)])
            
            drawnow
            pause(0.02);
            
        end
    end
end
end