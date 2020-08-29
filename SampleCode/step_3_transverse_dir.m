function step_3_transverse_dir
% This file computes the optimal transverse coordinates as defined the
% direction of the normal vector, z(\tau).

% set params
PFL = false;
coordSet = 'abs';
unactuated_nominal_trajectory = true;
saveFolder = 'alpha_3_abs_coord_unactuated';

if ~exist(saveFolder,'dir')
    mkdir(saveFolder)
end
 
if PFL
    switch coordSet
        case 'rel'
            loadFolder = 'alpha_2_2_closed_loop_PFL';
            loadedValues = load(sprintf('%s/results_alpha_2_2_PFL.mat',loadFolder));
        case 'abs'
            loadFolder = 'alpha_2_2_closed_loop_PFL_abs';
            loadedValues = load(sprintf('%s/results_alpha_2_2_PFL.mat',loadFolder));
    end
    nominal_PFL_u_pp = loadedValues.nominal_PFL_u_pp;
    f_PFL = loadedValues.f_PFL;
else
    % loadFolder = 'alpha_2_1_closed_loop_full';
    % loadFolder = 'Unactuated_alpha_2_1_closed_loop_full';
    switch coordSet
        case 'rel'
            loadFolder = 'alpha_2_1_closed_loop_full_rel_coord';
            if unactuated_nominal_trajectory
                loadedValues = load(sprintf('%s/result_0_2_closed_loop_sim_no_input.mat',loadFolder));
            else
                loadedValues = load(sprintf('%s/result_0_2_closed_loop_sim_f_full_feedback.mat',loadFolder));
            end
        case 'abs'
            if unactuated_nominal_trajectory
                loadFolder = 'alpha_2_1_closed_loop_full_abs_coord_unactuated';
                loadedValues = load(sprintf('%s/result_0_2_closed_loop_sim_no_input.mat',loadFolder));
            else
                loadFolder = 'alpha_2_1_closed_loop_full_abs_coord';
                loadedValues = load(sprintf('%s/result_0_2_closed_loop_sim_f_full_feedback.mat',loadFolder));
            end
    end
    feedbackLaw = loadedValues.feedbackLaw;

end

limCycleLength = loadedValues.limCycleLength;

nominal_u_pp = loadedValues.nominal_u_pp;

coordType = loadedValues.coordType;
controllerType = loadedValues.controllerType;
symImpactMap = loadedValues.symImpactMap;
foot_impact = loadedValues.foot_impact;

p = loadedValues.p;
f_full = loadedValues.f_full;
xDesiredpp = loadedValues.xDesiredpp;
ode45_max_time = loadedValues.ode45_max_time;
%% closed loop dynamics


if PFL
    switch coordType
        case 'rel'
            f_closed_loop = @(x) f_PFL_rel(x, input_PFL_relCoord_PFL_f(x, input_PD_FL_relCoord_full_f(x,p), p),p);
        case 'abs'
            f_closed_loop = @(x) f_PFL_abs(x, input_PFL_absCoord_PFL_f(x, input_PD_FL_abs_full_f(x,p), p),p);
    end
else
    f_closed_loop = @(x) f_full(x, feedbackLaw(x,p));
end


for normalDelta=[0.06]
for smoothCostCoef=[0]
for dz2dtCoef=[0]
for dotProdCoef=[0]
% for pOpt = [8,10,12,20,40,60,80]
for pOpt = [100]
    timer = tic();
% for normalDelta=[0.1]
% % for normalDelta=[0.3,0.4]
% for smoothCostCoef=[85]
% % for dz2dtCoef=[2,4,8]
% for dz2dtCoef=[48]
% for dotProdCoef=[0]
    
    % see if this combination has already been computed, if so, skip!
    try
        loadedVar = load(sprintf('zResults_2015_08_04_8norm_normalMin_%.2f_smoothCost_%d_dz2dtCoef_%d_dotProdCoef_%d.mat',normalDelta,smoothCostCoef,dz2dtCoef,dotProdCoef));
        display(sprintf('Already found zDir for normalMin=%.2f_smoothCost_%d_dz2dtCoef_%d_dotProdCoef_%d.mat',normalDelta,smoothCostCoef,dz2dtCoef,dotProdCoef))
        clear loadedVar
        continue
        
    catch err
        display(sprintf('Computing zDir for normalMin=%.2f  smoothCost=%d  dz2dtCoef=%d  dotProdCoef=%d',normalDelta,smoothCostCoef,dz2dtCoef,dotProdCoef))
                    
    end
    
%% Parameters
bezOrder = 5;

%% State variable initial condition
% q = loadedValues.initial_conditions'; % swing; stance; swingDot; stanceDot
q = xDesiredpp(0);

%% Simulation settings
number_steps = 8;  	%number of steps of the compass gait to simulate.

%% Switching event

% Use matlab's built in numerical solver to calculate the
% trajectory of the continuous phase. The continuous phase stops when the
% next foot touches the ramp. Set the 'Events' option so that the
% solver can calculate the impact point to a high precision.

odeParam = odeset('Events',foot_impact, 'RelTol',4.3e-9); % termination condition/settings

%% Figure
% figure()
% title('Phase Portrait of compass gait')
% xlabel('\theta [rad]')
% ylabel('d\theta/dt [rad/s]')
% hold on

% create (1 x num_steps) struct array with fields 't'(time vector) and
% 'q'(trajectory vector) for each step.
Step = struct('t',cell(1,number_steps),'q',cell(1,number_steps));

%% Simulation
for i = 1:number_steps
        
    % State vector: q = [theta(1)  theta(2) theta_dot(1) theta_dot(2)]
    %
    % Each column of q is state and each row corresponds to a time step.
    % now, q' is 'initial condition'
    [Step(i).t,Step(i).q] = ode45(@(t,q) f_closed_loop(q), [0,ode45_max_time], q', odeParam);
    
    % Calculate alpha, the angle between the legs.
    alpha = abs((Step(i).q(end,1) - Step(i).q(end,2)))/2;
    
    % Update q to the state in the instant before impact.
    q_final = Step(i).q(end,:);
    
    % impact map
    [deltaq, deltaqDot] = symImpactMap(q_final,p);
    % Update q to the state in the instant after foot impact.
    q = [deltaq*[q_final(1:2)']; deltaqDot*[q_final(3:4)]']';
    
    % Plot trajectory of one leg.
    
    if i>0  % start plotting after the i-th step.
        
        if mod(i,2) == 0
            % plot continuous phase (stance)
%             plothandle1 = plot(Step(i).q(:,2),Step(i).q(:,4),'color','blue');
            % plot discrete phase
%             linehandle1 = line([Step(i).q(end,2) q(1)],[Step(i).q(end,4) q(3)],'LineStyle','- -','color','blue');
            
            % take 'stable' limit cycle
            if i>=(number_steps-2)
                % highlight 'stable stance leg' in red
%                 set(plothandle1,'color','red')
%                 set(linehandle1,'color','red')
                
                % record the limit cycle
                stableStance = Step(i);
            end
        else
            % plot continuous phase (swing)
%             plothandle2 = plot(Step(i).q(:,1),Step(i).q(:,3),'black');
            % plot discrete phase
%             linehandle2 = line([Step(i).q(end,1) q(2)],[Step(i).q(end,3) q(4)],'LineStyle','- -','color','black');
            
            % take 'stable' limit cycle
            if i>=(number_steps-2)
                % highlight 'stable swing leg' in green
%                 set(plothandle2,'color','green')
%                 set(linehandle2,'color','green')
                
                
            end
        end
    end
end

%% Sample points along trajectories

stateSize = length(q);

sampleNo = 800;
decisionVarSize = (stateSize/2)*sampleNo;  % only use geometric values here...
% normalDelta = 0.4;
% smoothCostCoef = 10;
stanceTimeVec = stableStance.t;
stancePeriod = stanceTimeVec(end) - stanceTimeVec(1);
stance0 = stableStance.q(1,:);

[sampleStance.t,sampleStance.q] = ode45(@(t,q) f_closed_loop(q), linspace(0,stancePeriod,sampleNo), stance0', odeParam);
sampleTimeSpacing = stancePeriod/(sampleNo-1);
sampleTimeVec = sampleStance.t';

% plot(sampleStance.q(:,2),sampleStance.q(:,4),'x','color','red')

% plot swing leg samples
% plot(sampleStance.q(:,1),sampleStance.q(:,3),'x')

%% Transverse
% z direction

% sample f along sample points
fSample = zeros(4,sampleNo);
fSampleNormalised = zeros(4,sampleNo);
for iter = 1:sampleNo
    % compute f
    fSample(:,iter)=f_closed_loop(sampleStance.q(iter,:));
    % normalise f
    fSampleNormalised(:,iter) = fSample(:,iter)/norm(fSample(:,iter));
end

switch coordType
    case 'rel'
        impactSurfaceZ_start = [-1/sqrt(3) 2/sqrt(3) 0 0];
        impactSurfaceZ_end = [-1/sqrt(3) 2/sqrt(3) 0 0];
    case 'abs'
        impactSurfaceZ_start = [ -1/sqrt(2), -1/sqrt(2), 0 ,0];
        impactSurfaceZ_end = [ -1/sqrt(2), -1/sqrt(2), 0 ,0];
end

PStart = atan2(impactSurfaceZ_start(2), impactSurfaceZ_start(1));
PEnd = atan2(impactSurfaceZ_end(2), impactSurfaceZ_end(1));

% extract the f at control points as initial seeds
bezOrder = 5;

controlPtIndex = floor(linspace(1,sampleNo,bezOrder+1)); % corresponding index between control pts for Bezier and fSample
P0Index = controlPtIndex(2:end-1); % get rid of first and last point to obtain the index for initial seeds of P0
fAtP0Index = fSample(:,P0Index);
P0 = atan2(fAtP0Index(2,:),fAtP0Index(1,:));

zDecisionSize = stateSize/2;

% set up optimisation of transversal surfaces
lb = [];
ub = [];

%objective(z0)
warning('off','all')
optOptions = optimoptions('fmincon','Display','final-detailed'); % default MaxFunEvals:16000  MaxIter:400
[optP,fval,exitflag,optOutput] = fmincon(@objective,P0,     [],    [],      [],    []   ,lb,ub,@optEqCon,optOptions);
zMatrix= zPolyThetaAtTau(optP,sampleTimeVec);
zDotMatrix = gradient(zMatrix,sampleTimeSpacing);
optP;
finished_time = toc(timer);

%% compute piecewise polynomial Pi, PiDot
zSpline = spline(sampleStance.t,zMatrix);
zDotSpline = fnder(zSpline);
xStarSpline = spline(sampleStance.t,sampleStance.q');

bezPtsAll = [PStart, optP, PEnd];
xStarpp = @(t) ppval(xStarSpline,t);
zpp = @(t) [cos(bezier(bezPtsAll,t/stancePeriod)); sin(bezier(bezPtsAll,t/stancePeriod)); 0; 0];
% zDotpp = @(t) ppval(zDotSpline,t);
zDotpp = @(t) [-sin(bezier(bezPtsAll,t/stancePeriod))*bezierDer1(bezPtsAll,t/stancePeriod); cos(bezier(bezPtsAll,t/stancePeriod))*bezierDer1(bezPtsAll,t/stancePeriod); 0; 0];

Pipp = @(t) [-sin(bezier(bezPtsAll,t/stancePeriod)), cos(bezier(bezPtsAll,t/stancePeriod)), 0, 0; 0 0 1 0; 0 0 0 1]; 
PiDotpp = @(t) [-cos(bezier(bezPtsAll,t/stancePeriod))*bezierDer1(bezPtsAll,t/stancePeriod), -sin(bezier(bezPtsAll,t/stancePeriod))*bezierDer1(bezPtsAll,t/stancePeriod), 0, 0; 0 0 0 0; 0 0 0 0]; 

clear loadedVar
eval(sprintf('save %s/zResults_2015_08_04_8norm_normalMin_%.2f_smoothCost_%d_dz2dtCoef_%d_dotProdCoef_%d_pOpt_%d.mat',saveFolder,normalDelta,smoothCostCoef,dz2dtCoef,dotProdCoef,pOpt))

end
end
end
end
end

% optimisation equality and inequality constratints (norm=1 constraint)
    function [c, ceq] = optEqCon(PIn) % requires c<=0, and ceq==0
        
        ceq = [];
        %% REQUIRE dot(z,f) to be greater than minNormal at all sampled pts
        % calculate z, given tau
        
        zMat = zPolyThetaAtTau(PIn,sampleTimeVec);
        
        zDotf = dot(fSampleNormalised,zMat);
        
        c=-zDotf+normalDelta; % need c<=0, so that zdotf>normalDelta
        
    end
    
% computes the bezier polynomial representation of theta at tau
    function zAtTau = zPolyThetaAtTau(P,tau)
        % INPUTS:
        % PVec: vector containing the control points for bezier polynomial defining theta
        % tau: vector containing tau at which z is to be evaluated
        % OUTPUTS:
        % returns zAtTau, where each column is the output of z at a given tau input
        period = sampleTimeVec(end);
        
        t = tau./period;
        
        tauSize=length(tau);
        
        oneMinust = 1-t;
        oneMinust5 = oneMinust.^5; oneMinust4=oneMinust.^4; oneMinust3=oneMinust.^3; oneMinust2=oneMinust.^2;
        
        thetaAtTau = (oneMinust5).*PStart + 5.*t.*oneMinust4.*P(1) + 10.*t.^2.*oneMinust3.*P(2) +...
            10.*t.^3.*oneMinust2.*P(3) + 5.*t.^4.*oneMinust.*P(4) + t.^5.*PEnd;
             
        zGeometricAtTau = [cos(thetaAtTau); sin(thetaAtTau)];
        
        zAtTau = [zGeometricAtTau; zeros(zDecisionSize, tauSize)];
    end
%% optimisation objective
    function obj = objective(PIn)
        
        zMat = zPolyThetaAtTau(PIn,sampleTimeVec);
        
        dzdt = gradient(zMat,sampleTimeSpacing);
        
%         dz2dt = diff(dzdt,2,2)./sampleTimeSpacing;
%         dz2dtDiff = diff(zMat,2,2);
%         dzdtDiff = diff(zMat,1,2);
        
        dzdtpnorm = sum(dzdt.^pOpt, 1);
        
%         innerpnorm = dot(zMat,fSample).^pOpt;
        innerpnorm = (sqrt(dot(zMat,fSample).*dot(zMat,fSample)).*dot(zMat,fSample)).^pOpt;
        
        integrand = sum(dzdtpnorm./innerpnorm)^(1/pOpt);
        obj = integrand;
    end
    
    
function sum = bezier(a,s)
        % Compute bezier polynomial output and derivative...
        % Compute the Bezier polynomial
        sum = 0;
        der1 = 0;
        der2 = 0;
        order = length(a) - 1;
        % compute bezier
        for j = 1:(order+1)
            k = j-1;
            sum = sum + a(j)*(s^k)*((1-s)^(order-k))*factorial(order)/(factorial(k)*factorial(order-k));
        end
    end

    function der1 = bezierDer1(a,s)
        % Compute bezier polynomial output and derivative...
        sum = 0;
        der1 = 0;
        der2 = 0;
        order = length(a) - 1;
        for j = 1:order
            k = j-1;
            der1 = der1 + (a(j+1) - a(j))*(factorial(order)/(factorial(k)*factorial(order-k-1)))*s^k*(1-s)^(order-k-1);
        end
    end    


end