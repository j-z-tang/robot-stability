function step_4_compute_lyapunov


global stableEpsilon QFull

QFullList={diag([100 1 10 0.1])};
% The user may try other Q values, such as the following
% QFullList={diag([100 1 10 0.1]),diag([1 100 0.1 10]), diag([1 1 1 1]), diag([50 1 5 0.5])};
loadFolder = 'alpha_3_abs_coord_unactuated';
saveFolderRiccati = 'alpha_4_lyap_abs_coord_unactuated';

if ~exist(saveFolderRiccati,'dir')
    mkdir(saveFolderRiccati)
end

for normalDeltaNow=[0.06]
for smoothCostCoefNow=[0]
for dz2dtCoefNow=[0]
for dotProdCoefNow=[0]
for stableEpsilon = [1];
    for discreteQFactor = [1]
    for QFullIndex=1:length(QFullList)
       pOpt = 100; 
        LoadedVar = load(sprintf('%s/zResults_2015_08_04_8norm_normalMin_%.2f_smoothCost_%d_dz2dtCoef_%d_dotProdCoef_%d_pOpt_%d.mat',loadFolder,normalDeltaNow,smoothCostCoefNow,dz2dtCoefNow,dotProdCoefNow,pOpt));
p = LoadedVar.p;
bezPtsAll = LoadedVar.bezPtsAll;
PFL = LoadedVar.PFL;
coordType = LoadedVar.coordType;
controllerType = LoadedVar.controllerType;
symImpactMap = LoadedVar.symImpactMap;
foot_impact = LoadedVar.foot_impact;

if PFL
    switch coordType
        case 'rel'
            nominal_u = LoadedVar.nominal_PFL_u_pp;
            dfdxNowpp = @(x1,x2,x3,x4,  u) dfdx_PFL_rel([x1,x2,x3,x4],u,p);
            dfduNowpp = @(x1,x2,x3,x4)  dfdu_PFL_rel([x1,x2,x3,x4],p);
            
            f = @(x, u)  f_PFL_rel(x,u,p);
        case 'abs'
            nominal_u = LoadedVar.nominal_PFL_u_pp;
            dfdxNowpp = @(x1,x2,x3,x4,  u) dfdx_PFL_abs([x1,x2,x3,x4],u,p);
            dfduNowpp = @(x1,x2,x3,x4)  dfdu_PFL_abs([x1,x2,x3,x4],p);
            
            f = @(x, u)  f_PFL_abs(x,u,p);
    end
    
else % NO PFL, full dynamics
    switch coordType
        case 'rel'
            nominal_u = LoadedVar.nominal_u_pp
            f_full = LoadedVar.f_full;
            dfdxNowpp = @(x1,x2,x3,x4, u) dfdx_full_rel([x1,x2,x3,x4],u,p);
            dfduNowpp = @(x1,x2,x3,x4)  dfdu_full_rel([x1,x2,x3,x4],p);
            
            f = f_full;
            
        case 'abs'
            nominal_u = LoadedVar.nominal_u_pp;
            f_full = LoadedVar.f_full;
            dfdxNowpp = @(x1,x2,x3,x4, u) dfdx_full_abs([x1,x2,x3,x4],u,p);
            dfduNowpp = @(x1,x2,x3,x4)  dfdu_full_abs([x1,x2,x3,x4],p);
            
            f = f_full;
    end
end



checkInlineFunc = false;

QFull=QFullList{QFullIndex};

if exist(sprintf('%s/z_Pi_P_bezierTheta_extra_stable_ep_%d__Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt%d_Qd_factor_%d.mat',saveFolderRiccati,stableEpsilon,mat2str(diag(QFull)'),normalDeltaNow,pOpt,discreteQFactor),'file')
    prompt = 'Previous computation found.  Recompute Riccati solution?  Note if the .mat files were moved from another machine, you must recompute them.  Y/N [N]';
    userInput = input(prompt,'s');
    if strcmp(userInput,'Y')
        display('Recomputing Riccati solution')
    else
        return 
    end
end

zMatrix = LoadedVar.zMatrix;
zDotMatrix = LoadedVar.zDotMatrix;

optP = LoadedVar.optP; % optimised bezier coefficients; not to be confused with pOpt (p-norm for optimization) above.
PStart = LoadedVar.PStart;
PEnd = LoadedVar.PEnd;
stancePeriod = LoadedVar.stancePeriod;
sampleTimeVec = LoadedVar.sampleTimeVec;
sampleStance = LoadedVar.sampleStance;

normalDelta = LoadedVar.normalDelta;
smoothCostCoef = LoadedVar.smoothCostCoef;
dz2dtCoef=LoadedVar.dz2dtCoef;
dotProdCoef=LoadedVar.dotProdCoef;

xStarpp = LoadedVar.xStarpp;

zpp = LoadedVar.zpp;
zDotpp = LoadedVar.zDotpp;

Pipp = LoadedVar.Pipp;
PiDotpp = LoadedVar.PiDotpp;

%% check piecewise polys against loaded values
if checkInlineFunc
    figure(); hold on; title('z(1) values')
    plot(sampleTimeVec, zMatrix(1,:),'o');
    zFrompp = zeros(4,length(sampleTimeVec));
    for index = 1:length(sampleTimeVec)
        zFrompp(:,index) = zpp(sampleTimeVec(index));
    end
    plot(sampleTimeVec, zFrompp(1,:));
    
    figure(); hold on; title('zDot(1) values')
    plot(sampleTimeVec, zDotMatrix(1,:));
    zDotFrompp = zeros(4,length(sampleTimeVec));
    for index = 1:length(sampleTimeVec)
        zDotFrompp(:,index) = zDotpp(sampleTimeVec(index));
    end
    plot(sampleTimeVec, zDotFrompp(1,:));
    
    [PippOld,PiDotppOld] = piTau(@(t)ppval(zSpline,t), @(t)ppval(zDotSpline,t)); % legacy code...
    figure(); hold on; title('PiDot(1,1) values')
    PiDot11Frompp = zeros(1,length(sampleTimeVec));
    PiDot11FromppNew = zeros(1,length(sampleTimeVec));
    for index = 1:length(sampleTimeVec)
        PiDotppNow = PiDotpp(sampleTimeVec(index));
        PiDot11Frompp(:,index) = PiDotppNow(1,1);
        PiDotppNewNow = PiDotppOld(sampleTimeVec(index));
        PiDot11FromppNew(:,index) = PiDotppNewNow(1,1);
    end
    plot(sampleTimeVec, PiDot11FromppNew);
    plot(sampleTimeVec, PiDot11Frompp)
end

%% compute Ricatti equation
% QFull = diag([100 50 10 5]);
noSamplePt = length(zMatrix);
xspanLimCycSample = sampleStance.q'; % each x is a column
tspanLimCycSample = sampleStance.t; % column 
limCycDuration = tspanLimCycSample(end);

% transformation matrix before and after
PiStart = Pipp(0);        PiEnd = Pipp(limCycDuration);
PiStartDot = PiDotpp(0);  PiEndDot = PiDotpp(limCycDuration);

% impact point
impactPt = xspanLimCycSample(:,end);
xSym = sym('xSym',[4,1]);
xSym = sym(xSym,'real');
qSw = xSym(1);
qSt = xSym(2);
qSwDot = xSym(3);
qStDot = xSym(4);

[delta_q, delta_qDot] = symImpactMap([qSw,qSt],p);

impactTransformMatrix = [delta_q, zeros(2);
                        zeros(2), delta_qDot];
                      
discreteJump = impactTransformMatrix * xSym;
discreteJumpdx = jacobian(discreteJump,xSym);
% transverse linearization of impact map
transLinImpactMap = PiStart * discreteJumpdx * PiEnd';
transLinImpactMapAtImpactPt = double(subs(transLinImpactMap, xSym, impactPt));

%% Solve the new P via the jump-Riccati differential equation
Pdimension = 3;
no_iter = 20; %number of times we should solve the ricatti equation to obtain periodic P
clear PIter
PIter = struct('t',cell(1,no_iter), 'PColumn',cell(1,no_iter));
% set initial condition
currentPMatrix = eye(Pdimension);
ep = 1e-2; % original

for i=1:no_iter
    display(sprintf('Computing solution to the Riccati equation - iteration no.%d',i))
    % jump ricatti condition
    %  display('jump condition:')
    currentPMatrix = transLinImpactMapAtImpactPt'*currentPMatrix*transLinImpactMapAtImpactPt +  PiEnd * ((1/discreteQFactor)*QFull) * PiEnd';
    % reshape P to a column for ode45
    currentPColumn = reshape(currentPMatrix,[],1);
    % solve continuous P ode45
   [PIter(i).t, PRowOutput] = ode45(@Ricatti,[limCycDuration 0],currentPColumn);
    PColumnOutput = PRowOutput';
    PIter(i).PColumn = PColumnOutput;
    %display('Final P value for this iteration:')
    currentPMatrix = reshape(PColumnOutput(:,end),Pdimension,Pdimension) ;   
    % save currentP, start again
end

% save P as spline
PColumnSpline = spline(PIter(i).t,PIter(i).PColumn);
PDotColumnSpline_fnder = fnder(PColumnSpline); % this is the old method, we found to be not so accurate...
PMatrixpp = @(t) reshape(ppval(PColumnSpline,t),Pdimension,Pdimension);
PDotMatrixpp_fnder = @(t) reshape(ppval(PDotColumnSpline_fnder,t),Pdimension,Pdimension);  % this is the old method (with fnder)

clear PDotColumnSpline_Ricatti

for iter = 1:length(PIter(i).t)
    PDotColumnSpline_Ricatti(:,iter) = Ricatti(PIter(i).t(iter),PIter(i).PColumn(:,iter)); % compute PDot directly from ricatti eq.
end
PDotColumnSpline_Ricatti_Spline = spline(PIter(i).t,PDotColumnSpline_Ricatti);
PDotMatrixpp = @(t) reshape(ppval(PDotColumnSpline_Ricatti_Spline,t),Pdimension,Pdimension); % do the inline func with ricatti equation...


clear LoadedVar;
eval(sprintf('save %s/z_Pi_P_bezierTheta_extra_stable_ep_%d__Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt%d_Qd_factor_%d.mat',saveFolderRiccati,stableEpsilon,mat2str(diag(QFull)'),normalDelta,pOpt,discreteQFactor))

    end
    end
end
end
end
end
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
        % Compute the Bezier polynomial
        sum = 0;
        der1 = 0;
        der2 = 0;
        order = length(a) - 1;
        % compute first derivative
        for j = 1:order
            k = j-1;
            der1 = der1 + (a(j+1) - a(j))*(factorial(order)/(factorial(k)*factorial(order-k-1)))*s^k*(1-s)^(order-k-1);
        end
    end    

    function der2 = bezierDer2(a,s)
        % Compute bezier polynomial output and derivative...
        % Compute the Bezier polynomial
        sum = 0;
        der1 = 0;
        der2 = 0;
        order = length(a) - 1;
        % compute second derivative
        for j = 1:order-1
            k = j-1;
            der2 = der2 + (a(j+2) - 2*a(j+1) + a(j))*(factorial(order)/(factorial(k)*factorial(order-k-2)))*s^k*(1-s)^(order-k-2);
        end
    end
%% function to compute transverse ricatti equation
    function pDoutColumn = Ricatti(t,PIn)

        % reshape P into matrix form
        PMatrix = reshape(PIn,Pdimension,Pdimension);
        
        xStarEval = xStarpp(t);
        uStarEval = nominal_u(t);
        
        zNow = zpp(t); zdotNow=zDotpp(t);
        fNow = f(xStarEval, uStarEval);
%         fNow = f(xStarEval(1),xStarEval(2),xStarEval(3),xStarEval(4), uStarEval);
        PiNow = Pipp(t); PiDotNow = PiDotpp(t);
%         dfdxNow = dfdxNowpp(xStarEval, uStarEval);
        dfdxNow = dfdxNowpp(xStarEval(1),xStarEval(2),xStarEval(3),xStarEval(4), uStarEval);
        
%         dfduNow = dfduNowpp(xStarEval);
        dfduNow = dfduNowpp(xStarEval(1),xStarEval(2),xStarEval(3),xStarEval(4));
        
        
        dtaudx = (zNow'*dfdxNow*PiNow' + zdotNow'*PiNow')/(zNow'*fNow);
        
        transLinA = (PiDotNow*PiNow' + PiNow*dfdxNow*PiNow' - PiNow*fNow*dtaudx) + eye(3)*stableEpsilon;
        
        dtaudu = zNow'*dfduNow/(zNow'*fNow);
        transLinB = PiNow*dfduNow - PiNow * fNow * dtaudu;
        
        % compute Pdot
        Qtransverse = Pipp(t) * QFull * Pipp(t)';
        RInverse = 100;
        Poutput = (-1)*(transLinA'*PMatrix + PMatrix*transLinA - PMatrix*transLinB*RInverse*transLinB'*PMatrix + Qtransverse);  % jump riccati
%        Poutput = (-1)*(transLinA'*PMatrix + PMatrix*transLinA + Qtransverse); % jump lyapunov
        pDoutColumn = reshape(Poutput, [],1);
        

        
    end
end