function step_6_verification
% This function performs the SOS verification to search for Lyapunov
% functions which guarantee orbital stability of the limit cycle.

QFullListSize = 1;
loadFolder_Taylor = 'alpha_5_taylor_dynamics_abs_unactuated';
for normalDeltaNow=[0.06]
    for stableEpsilonNow = [1];
        for discreteQFactor = [1]
            for QFullIndex=1:QFullListSize
                QFullList={diag([100 1 10 0.1])};
                %QFullList={diag([100 1 10 0.1]), diag([50 1 5 0.5]), diag([1 1 1 1]) ,diag([1 100 0.1 10]) };  % can also try these
                if QFullIndex > length(QFullList)
                    continue;
                end
                QFullNow=QFullList{QFullIndex};
                pOpt = 100;
                %% Make folder for saving results
                loaded = load(sprintf('%s/spline_full_dynamics_taylor_stableEp_%d_Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d.mat',loadFolder_Taylor,stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor));
                if ~isequal(stableEpsilonNow,loaded.stableEpsilonNow) || ~isequal(normalDeltaNow,loaded.normalDeltaNow)
                    error('file corrupted')
                end
                QFull = loaded.QFull;
                pOpt = loaded.pOpt;
                saveFolder = sprintf('alpha_6_ABS_COORD_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
                
                if ~exist(saveFolder,'dir')
                    mkdir(saveFolder)
                end
                
                %% Get taylor series of transverse dynamics at this point into msspoly
                thisSectionFile = sprintf('%s/results_5_2_PRELIM_new_zafter_taylor_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d.mat',saveFolder,stableEpsilonNow,mat2str(diag(QFull)'),normalDeltaNow,pOpt);
                runTaylorSection = true;
                if exist(thisSectionFile,'file')
                    prompt = 'Previous computation found.  Recompute conversion to msspoly?  Note if the .mat files were moved from another machine, you must recompute them.  Y/N [N]';
                    userInput = input(prompt,'s');
                    if strcmp(userInput,'Y')
                        display('Recomputing Taylor expansions...')
                    else
                        runTaylorSection = false;
                    end
                end
                
                if runTaylorSection
                    RInverse = loaded.RInverse;
                    sampleNo = loaded.sampleNo;
                    sampleTimeVec = loaded.sampleTimeVec;
                    PMatrixpp = loaded.PMatrixpp;
                    PDotMatrixpp = loaded.PDotMatrixpp;
                    Pipp = loaded.Pipp;
                    PiDotpp = loaded.PiDotpp;
                    zpp = loaded.zpp;
                    zDotpp = loaded.zDotpp;
                    xStarpp = loaded.xStarpp;
                    f = loaded.f;
                    fFullTaylorSample = loaded.fFullTaylorSample;
                    dfduNowpp = loaded.dfduNowpp;
                    
                    limCycDuration = loaded.limCycDuration;
                    
                    nominal_u = loaded.nominal_u;
                    deltaq = loaded.deltaq;
                    deltaqDotTaylor = loaded.deltaqDotTaylor;
                    %% setup variables
                    display('The time now is')
                    datestr(clock,'yyyy-mm-dd_HH-MM-SS')
                    mosekopt
                    
                    totalSamples = 40;
                    N=totalSamples-1;
                    info = '';
                    display('=========================================================================================================')
                    %     display(['rho = ', num2str(rho)])
                    
                    tstart = tic;
                    
                    
                    L1 = cell([totalSamples 1]); L2 = cell([totalSamples 1]);
                    fTaylorAll = cell([totalSamples 1]);
                    solAll = cell([totalSamples 1]);
                    precomp = cell([totalSamples 1]);
                    
                    clear prog
                    %     prog = spotsosprog; % set up program, add variables to it later
                    xper = msspoly('xper',3); variables.xper = xper; % [prog,xper]=prog.newIndeterminate('xper',3);
                    x = msspoly('x',4); variables.x = x;  %[prog,x]=prog.newIndeterminate('x',4);
                    uCon = msspoly('ucon',1); variables.uCon = uCon; % [prog,uCon]=prog.newIndeterminate('ucon',1);
                    tau = msspoly('tau',1); variables.tau = tau; %  [prog,tau]=prog.newIndeterminate('tau',1);
                    xper1=xper(1); xper2=xper(2); xper3=xper(3); % for use when converting taylor series from symbolic maths to msspoly
                    x1=x(1); x2=x(2); x3=x(3); x4=x(4);
                    fTaylor = msspoly('ftay',4); %[prog,fTaylor]=prog.newFree(4);
                    taudot = msspoly('taud',1); %[prog,taudot] = prog.newFree(1);
                    g = msspoly('g',3); %[prog,g]=prog.newFree(3);
                    
                    fTaylor=fTaylor-fTaylor; g=g-g; taudot=taudot-taudot;
                    
                    % switch conditions
                    variables.deltaq = deltaq;
                    noElementsDeltaDot = numel(deltaqDotTaylor);
                    deltaDotColumn = reshape(deltaqDotTaylor,noElementsDeltaDot,1);
                    
                    deltaqDotColumnSpot = msspoly('delt',noElementsDeltaDot); deltaqDotColumnSpot=deltaqDotColumnSpot-deltaqDotColumnSpot;
                    for elementIndex=1:length(deltaDotColumn)
                        [elementCoef, elementMono] = coeffs(deltaDotColumn(elementIndex));
                        elementCoefDouble = double(elementCoef);
                        % iterate through each monomial
                        for monoIndex = 1:length(elementMono)
                            % fperp(i) = fperp(i) + char(monomials(monoIndex))*coeffs
                            %eval(['fperp(elementIndex)=fperp(elementIndex)+(',char(elementMono(monoIndex)),'*(',char(elementCoef(monoIndex)),'));'])
                            eval(['deltaqDotColumnSpot(elementIndex)=deltaqDotColumnSpot(elementIndex)+(',char(elementMono(monoIndex)),'*elementCoefDouble(monoIndex));'])
                        end
                    end
                    variables.deltaqDotColumnSpot = deltaqDotColumnSpot;
                    
                    display('Converting Taylor series to Spotless format...')
                    for sampleIndex=1:totalSamples
                        fTaylor=fTaylor-fTaylor; % reset fTaylor to 0
                        for elementIndex=1:length(fFullTaylorSample{sampleIndex})
                            [elementCoef, elementMono] = coeffs(fFullTaylorSample{sampleIndex}(elementIndex));
                            elementCoefDouble = double(elementCoef);
                            % iterate through each monomial
                            for monoIndex = 1:length(elementMono)
                                % fperp(i) = fperp(i) + char(monomials(monoIndex))*coeffs
                                %eval(['fperp(elementIndex)=fperp(elementIndex)+(',char(elementMono(monoIndex)),'*(',char(elementCoef(monoIndex)),'));'])
                                eval(['fTaylor(elementIndex)=fTaylor(elementIndex)+(',char(elementMono(monoIndex)),'*elementCoefDouble(monoIndex));'])
                            end
                        end
                        fTaylorAll{sampleIndex} = fTaylor;
                    end
                    clear fTaylor sampleIndex loaded;
                    save(thisSectionFile)
                else
                    load(thisSectionFile)
                end
                
                
                %% Compute preliminaries
                thisSectionFile = sprintf('%s/results_5_2_PRELIM_after_precomp_taylor__stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d.mat',saveFolder,stableEpsilonNow,mat2str(diag(QFull)'),normalDeltaNow,pOpt);
                
                runPrelimSection = true;
                if exist(thisSectionFile,'file')
                    prompt = 'Previous computation found.  Recompute preliminaries?  Note if the .mat files were moved from another machine, you must recompute them.  Y/N [N]';
                    userInput = input(prompt,'s');
                    if strcmp(userInput,'Y')
                        display('Recomputing Taylor expansions...')
                    else
                        runPrelimSection = false;
                    end
                end
                
                if runPrelimSection
                    ep = 1e-6;
                    
                    wDeg = 6; % degrees of polynomial in W
                    pdeg = 6; % degrees of polynomials for the multipliers
                    PDummyVec = msspoly('pdum',6); %[prog, PDummyVec] = prog.newIndeterminate('pdum',6);
                    PDummy = mss_v2s(PDummyVec);
                    V=xper.'*PDummy*xper;
                    
                    % dVdx
                    dVdxper = diff(V,xper);
                    % dVdP
                    dVdP = diff(V,PDummyVec); % row of [dVdP11, dVdP12.... dVdP33]
                    % metric monomials
                    metricMono = monomials([xper;tau],0:wDeg);
                    % matrix multiplier monomials
                    multiMono = monomials(xper,0:pdeg);
                    % s-procedure monomials
                    sMono = monomials(xper,0:pdeg);
                    dts=diff(sampleTimeVec);
                    % start precomp to compute z, zDot, f, Pi, PiDot, P PDot
                    display('computing z, P, Pi and its derivatives')
                    for i = 1:totalSamples
                        precomp{i}=struct();
                        timeNow = sampleTimeVec(i); precomp{i}.timeNow = timeNow;
                        
                        xStarNow = xStarpp(timeNow); precomp{i}.xStarNow = xStarNow;
                        uStarNow = nominal_u(timeNow); precomp{i}.uStarNow = uStarNow;
                        fNow = f(xStarNow, uStarNow); precomp{i}.fNow = fNow;
                        %         fNow = f(xStarNow(1),xStarNow(2),xStarNow(3),xStarNow(4), uStarNow); precomp{i}.fNow = fNow;
                        PiNow = Pipp(timeNow); precomp{i}.PiNow = PiNow;
                        PiDotNow = PiDotpp(timeNow); precomp{i}.PiDotNow = PiDotNow;
                        zNow = zpp(timeNow); precomp{i}.zNow = zNow;
                        zDotNow = zDotpp(timeNow); precomp{i}.zDotNow = zDotNow;
                        
                        
                        % get P from ricatti
                        PMatrixNow = PMatrixpp(timeNow); precomp{i}.PMatrixNow = PMatrixNow;
                        
                        PSymVec = mss_s2v(PMatrixNow); % symmetric part only
                        
                        PdotMatrixNow = PDotMatrixpp(timeNow); precomp{i}.PdotMatrixNow = PdotMatrixNow;
                        
                        PdotSymVec = mss_s2v(PdotMatrixNow); % symmetric part only
                        dVdtauNow = dVdP*PdotSymVec;  % dVdP = [dVdP11, dVdP12.... dVdP33]
                        precomp{i}.dVdtauNow = dVdtauNow;
                        
                        dVdxperNow = subs(dVdxper,PDummyVec,PSymVec); precomp{i}.dVdxperNow = dVdxperNow;
                        VNow = subs(V, PDummyVec,PSymVec); precomp{i}.VNow = VNow;
                        precomp{i}.RInverse = RInverse;
                    end
                    clear i timeNow xStarNow fNow PiNow PiDotNow zNow zDotNow PMatrixNow PdotMatrixNow dVdtauNow dVdxperNow VNow RInverse loaded
                    save(thisSectionFile)
                else
                    load(thisSectionFile)
                    clear loaded
                end
                
                
                
                %% Begin verification
                for SDSOS_verification=[1]
                    finishedFileName = sprintf('%s/result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_SDSOS_%d.mat',saveFolder,stableEpsilonNow,mat2str(diag(QFull)'),normalDeltaNow,pOpt,SDSOS_verification);
                    
                    runVerificationSec = true;
                    if exist(finishedFileName,'file')
                        prompt = 'Previous computations found.  Recompute SOS verification?  Note if the .mat files were moved from another machine, you must recompute them.  Y/N [N]';
                        userInput = input(prompt,'s');
                        if strcmp(userInput,'Y')
                            display('Recomputing Taylor expansions...')
                        else
                            runVerificationSec = false;
                        end
                    end
                    
                    if runVerificationSec %HAS NOT SUCCEEDED BEFORE, TRY NOW
                        diaryName = sprintf('%s/command_line_output_SDSOS_%d.txt',saveFolder,SDSOS_verification);
                        if exist(diaryName,'file')
                            delete(diaryName);
                        end
                        diary(diaryName)
                        diary on
                        display(sprintf('The time now is %s',datestr(clock,'yyyy-mm-dd_HH-MM-SS')))
                        
                        try
                            SDSOS_verification;
                            rhoTypeAll = {'poly_tau','piecewise_linear','chebyshev_poly'};
                            variables.rhoType = rhoTypeAll{3};
                            
                            variables.solver = @spot_mosek;
                            
                            variables.findMultiplierBalance = 1;
                            variables.optimizeRhoBalance = 1;
                            
                            if SDSOS_verification
                                variables.useSDSOS_multiplier_step = true; 
                                variables.useSDSOS_rho_step = true;
                            else
                                variables.useSDSOS_multiplier_step = false;
                                variables.useSDSOS_rho_step = false;
                            end
                            variables.discountRho_in_multiplier_step = true; % the maximised \rho in \rho step may not be strictly in the feasible set; discounting it helps numerically
                            variables.maxGamma_multiplier_step = false;
                            variables.multiplierOrder = 6;
                            
                            variables.periodicRho_rho_step = false;
                            rhoStepObjectiveAll = {'max_integral','max_single_slack','max_many_slack'};
                            variables.rhoStepObjective = rhoStepObjectiveAll{1};
                            
                            variables.cleanSOSPoly = true;
                            variables.SOSCleanTol = 1e-9;
                            
                            variables.negativity_buffer_epsilon_multiplier_step = 1e-6;
                            variables.negativity_buffer_DISCRETE_multiplier_step = 1e-6;
                            variables.negativity_buffer_epsilon_rho_step = 1e-5;
                            variables.negativity_buffer_DISCRETE_rho_step = 5e-4;
                            variables.plot_rho_step = false;
                            
                            variables.dfduNowpp = dfduNowpp;
                            
                            variables.limCycDuration = limCycDuration;
                            
                            
                            % fix rho, search multiplier
                            switch variables.rhoType
                                case 'poly_tau'
                                    rhoDeg = 20;   variables.rhoDeg=rhoDeg;
                                    variables.tauMono = monomials(tau,[0:rhoDeg]);
                                    rhoInit = 0.01;
                                    rhodotSol=[];
                                    rhoSol = [rhoInit;zeros(rhoDeg,1)];  % these are coefficients of rho
                                case 'piecewise_linear'
                                    rhoInit = 0.1;
                                    rhoSol = ones(40,1)*rhoInit;
                                    rhodotSol = zeros(39,1);
                                    variables.rhof = 0.1; % final rho
                                case 'chebyshev_poly'
                                    rhoDeg = 20;   variables.rhoDeg=rhoDeg;
                                    %rhoInit = 0.01
                                    rhoInit = 0.005;
                                    rhoSol = [rhoInit;zeros(rhoDeg,1)];
                                    rhodotSol=[];
                            end
                            
                            % fix multiplier, max rho
                            optimise_rho_mode = {'max_rho','rho_feasibility'};
                            
                            % fix rho, search multiplier
                            totalIteration = 6;
                            
                            % clear these variables from the previous case...
                            clear multiplierConds LSolStruct TAllStruct TiAllStruct rhoConds rhoSolStruct rhodotSolStruct rhointegralAll rho_feasibility_damping_all multiplier_step_rho_damping_all multiplier_step_time rho_step_time
                            for iter=1:totalIteration
                                display(sprintf('Starting iteration no. %d for rho-step and multipliers',iter))
                                
                                clear LSolNew
                                display('Beginning multiplier-step for this iteration...')
                                %         [LSolNew, MulsolAll, TAll, TiAll, SOSConds] = findMultiplier(precomp, variables, rhoSol, rhodotSol, fTaylorAll, multiMono, TAll, TiAll);
                                variables.discountRhoFactor = 1;
                                variables;
                                multiplier_feasible = false;
                                while (~multiplier_feasible)  &&  (variables.discountRhoFactor > 0)
                                    try
                                        display(sprintf('Running multiplier-step with a rho damping factor of %f',variables.discountRhoFactor));
                                        multiplier_step_clock = tic;
                                        [LSolNew, MulsolAll, TAll, TiAll, ~] = findMultiplier(precomp, variables, rhoSol, rhodotSol, fTaylorAll, multiMono); % balance every time
                                        multiplier_step_dur = toc(multiplier_step_clock);
                                        multiplier_feasible = true;
                                        multiplier_step_time(iter) = multiplier_step_dur
                                        multiplier_step_rho_damping_all(iter) = variables.discountRhoFactor;
                                        
                                    catch err
                                        display(sprintf('FAILED multiplier-step with damping factor = %f',variables.discountRhoFactor));
                                        variables.discountRhoFactor = variables.discountRhoFactor - 0.1;
                                        display(sprintf('Reducing multiplier-step damping factor to %f now', variables.discountRhoFactor))
                                    end
                                end
                                LSolStruct{iter} = LSolNew;
                                TAllStruct{iter} = TAll;
                                TiAllStruct{iter}=TiAll;
                                
                                display('Beginning rho-step for this iteration...')
                                % max rho first
                                [rhoSol, rhodotSol, rhointegral, ~] = optimizeRho(precomp, variables, LSolNew, fTaylorAll, dts, N, TAll, TiAll,optimise_rho_mode{1});
                                
                                % solve feasibility problems
                                variables.rho_feasibility_step_damping = 0.9; % discount factor
                                variables;
                                rho_feasible = false;
                                while ~rho_feasible &&  (variables.rho_feasibility_step_damping > 0)
                                    try
                                        display('To recover a strictly feasible rho, we now re-run rho-step as a feasibility problem.')
                                        display(sprintf('We set a lower-bound on the previously optimized objective that is dampled by a factor of %f ...', variables.rho_feasibility_step_damping))
                                        rho_step_clock = tic;
                                        [rhoSol, rhodotSol, rhointegral, ~] = optimizeRho(precomp, variables, LSolNew, fTaylorAll, dts, N, TAll, TiAll,optimise_rho_mode{2}, rhointegral);
                                        rho_step_dur = toc(rho_step_clock);
                                        rho_step_time(iter) = rho_step_dur
                                        rho_feasible = true;
                                        rho_feasibility_damping_all(iter) = variables.rho_feasibility_step_damping;
                                    catch err
                                        display(sprintf('MOSEK was unable to recover a strictly feasible rho from previously optimized objective with damping factor = %f... We now reduce the damping factor further...',variables.rho_feasibility_step_damping))
                                        variables.rho_feasibility_step_damping = variables.rho_feasibility_step_damping*0.8;
                                        % display(sprintf('Trying again to recover strictly feasible rho by reducing the optimized objective damping factor to %f now', variables.rho_feasibility_step_damping))
                                    end
                                end
                                rhoSol;
                                rhoSolStruct{iter} = rhoSol;
                                rhodotSolStruct{iter} = rhodotSol;
                                rhointegralAll(iter,1) = rhointegral
                                
                                if iter > 1
                                    if (abs(rhointegral - rhointegralAll(end-1))/rhointegralAll(end-1)) < 0.01 % assume converged if we are within 1 percent
                                        break
                                    end
                                end
                                % Save partial result:
                                % save(sprintf('%s/result_5_2_bilinear_PARTIAL_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_SDSOS_%d.mat',saveFolder,stableEpsilonNow,mat2str(diag(QFull)'),normalDeltaNow,pOpt,SDSOS_verification))
                            end
                            save(finishedFileName)
                            diary off
                            clear multiplierConds LSolStruct TAllStruct TiAllStruct rhoConds rhoSolStruct rhodotSolStruct rhointegralAll rho_feasibility_damping_all
                        catch err
                            display('SOMETHING WENT WRONG, MOVING ON TO NEXT EXAMPLE')
                            display(sprintf('The time now is %s',datestr(clock,'yyyy-mm-dd_HH-MM-SS')))
                            diary off
                        end
                        clear multiplierConds LSolStruct TAllStruct TiAllStruct rhoConds rhoSolStruct rhodotSolStruct rhointegralAll
                    end % end if ~exist file
                    
                end % end iterating SD/SDSOS
                
            end
        end
    end
end
end



function [LSolNew, solAll, TAll, TiAll, SOSConds] = findMultiplier(precomp, variables, rhoAll, rhodot, fTaylorAll, multiMono, TAll, TiAll)

xper=variables.xper;
x=variables.x;
uCon=variables.uCon;
tau=variables.tau;

display('Beginning iteration for finding multipliers')

SOSConds = cell(length(precomp),1);

if variables.discountRho_in_multiplier_step
    rhoAll = rhoAll*variables.discountRhoFactor;
end

%% switching condition
deltaqDotColumnSpot = variables.deltaqDotColumnSpot;
deltaqDot = reshape(deltaqDotColumnSpot,2,2);
deltaq = variables.deltaq;

prog = spotsosprog;
prog = prog.withIndeterminate(xper);

sampleIndex = 100;
switch variables.rhoType
    case 'poly_tau'
        %     prog = prog.withPos(subs(rhoPoly, tau, precomp{end}.timeNow) - variables.rhof);
        % need  V0(t+)/rho(t+) < V0(t-)/rho(t-)   i.e. V0(t-)rho(t+) - V0(t+)rho(t-) = SOS
        tauMono = variables.tauMono;
        rhoMinus = subs(rhoAll'*tauMono, tau, precomp{end}.timeNow);
        rhoPlus = subs(rhoAll'*tauMono, tau, precomp{1}.timeNow);
        
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
    case 'piecewise_linear'
        
        rhoMinus = rhoAll(end);
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        rhoPlus = rhoAll(1);
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
    case 'chebyshev_poly'
        
        rhoMinus = chebyshevT(0:variables.rhoDeg, precomp{end}.timeNow) * rhoAll;
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        rhoPlus = chebyshevT(0:variables.rhoDeg, precomp{1}.timeNow) * rhoAll;
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
end
switchingCondition = V0Minus*rhoPlus - V0Plus*rhoMinus ; % needs to be > 0
VMinus = precomp{end}.VNow/double(rhoMinus);

% balancing
if nargin < 7
    S1=.5*double(subs(diff(diff(VMinus,xper)',xper),xper,0*xper));
    S2=.5*double(subs(diff(diff(switchingCondition,xper)',xper),xper,0*xper));
    [T,D] = balanceQuadForm(S1,S2);
    TAll{sampleIndex}=T;
    Ti = inv(T);
    TiAll{sampleIndex}=Ti;
else
    display('use existing T and Ti for switching')
end

xperBal = TAll{sampleIndex}*xper;
switchingConditionBal = subs(switchingCondition, xper, xperBal);

% define multipliers to search over
[prog, LNew{sampleIndex}{1}] = newFreePoly(prog,monomials(xper,0:variables.multiplierOrder));

regionIneq = rhoMinus - xperBal'*(precomp{end}.PMatrixNow)*xperBal;
% [prog,gamma] = newPos(prog,1);
buffer = variables.negativity_buffer_DISCRETE_multiplier_step;

if variables.useSDSOS_multiplier_step
    if variables.cleanSOSPoly
        prog = prog.withSDSOS(clean(-buffer*(xper.'*xper) + switchingConditionBal - LNew{sampleIndex}{1}*regionIneq,variables.SOSCleanTol));
        prog = prog.withSDSOS(LNew{sampleIndex}{1});
    else
        prog = prog.withSDSOS(-buffer*(xper.'*xper) + switchingConditionBal - LNew{sampleIndex}{1}*regionIneq);
        prog = prog.withSDSOS(LNew{sampleIndex}{1});
    end
else
    if variables.cleanSOSPoly
        prog = prog.withSOS(clean(-buffer*(xper.'*xper) + switchingConditionBal - LNew{sampleIndex}{1}*regionIneq,variables.SOSCleanTol));
        prog = prog.withSOS(LNew{sampleIndex}{1});
    else
        prog = prog.withSOS(-buffer*(xper.'*xper) + switchingConditionBal - LNew{sampleIndex}{1}*regionIneq);
        prog = prog.withSOS(LNew{sampleIndex}{1});
    end
end
opt = spot_sdp_default_options();
opt.verbose = 0;

solSearchMultiplier = prog.minimize(msspoly(0),variables.solver,opt);
solAll{sampleIndex} = solSearchMultiplier;
if ~strcmp(char(solSearchMultiplier.info.status),'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
    display(sprintf('ISSUES FOR SAMPLE %d, %s', sampleIndex, solSearchMultiplier.info.status))
    
else
    
end

LSolNew{sampleIndex}{1} = subs(solSearchMultiplier.eval(LNew{sampleIndex}{1}),xper, TiAll{sampleIndex}*xper);



%% continuous conditions
switch variables.rhoType
    case 'poly_tau'
        iterLastSlide = 40;
    case 'piecewise_linear'
        
        %         iterLastSlide = 39;
        
    case 'chebyshev_poly'
        iterLastSlide = 40;
end

for sampleIndex = 1:iterLastSlide % last slide is switching
    clear prog % each multiplier is a new prog!
    display(sprintf('Computing multipliers for surface no. %d',sampleIndex))
    timeNow = precomp{sampleIndex}.timeNow;
    
    
    switch variables.rhoType
        case 'poly_tau'
            tauMono = variables.tauMono;
            rhoNow = double(subs(rhoAll'*tauMono, tau, timeNow));
        case 'piecewise_linear'
            
            rhoNow = rhoAll(sampleIndex);
            
        case 'chebyshev_poly'
            rhoNow = chebyshevT(0:variables.rhoDeg, timeNow) * rhoAll;
    end
    
    clear prog
    prog = spotsosprog;
    prog = prog.withIndeterminate(xper);
    
    VNow = precomp{sampleIndex}.VNow/rhoNow;
    
    [DV, taudotDenominator]=transverseConditions(sampleIndex,precomp{sampleIndex},rhoAll,rhodot,fTaylorAll, variables);
    
    % balancing
    if nargin < 7
        S1=.5*double(subs(diff(diff(VNow,xper)',xper),xper,0*xper));
        S2=.5*double(subs(diff(diff(DV,xper)',xper),xper,0*xper));
        [T,D] = balanceQuadForm(S1,S2);
        TAll{sampleIndex}=T;
        Ti = inv(T);
        TiAll{sampleIndex}=Ti;
    else
        % use existing TAll and TiAll passed in
        display('use existing T and Ti...')
    end
    
    if variables.findMultiplierBalance % balance!
        xperBal = TAll{sampleIndex}*xper;
    else
        xperBal = xper;
    end
    DVBal = subs(DV,xper,xperBal);
    VNowBal = subs(VNow,xper,xperBal);
    taudotDenominatorBal = subs(taudotDenominator, xper, xperBal);
    
    regionIneq = rhoNow - xperBal'*(precomp{sampleIndex}.PMatrixNow)*xperBal;
    
    % define multipliers to search over
    [prog, LNew{sampleIndex}{1}] = newFreePoly(prog,monomials(xper,0:variables.multiplierOrder));
    
    [prog, LNew{sampleIndex}{2}] = newFreePoly(prog,monomials(xper,0:variables.multiplierOrder));
    
    if variables.maxGamma_multiplier_step
        [prog,gamma] = newPos(prog,2);  SOSConds{sampleIndex}.gamma = gamma;
        gamma = gamma + variables.negativity_buffer_epsilon_multiplier_step;
    else
        gamma = zeros(2,1);
        gamma = gamma + variables.negativity_buffer_epsilon_multiplier_step;
    end
    
    if variables.useSDSOS_multiplier_step
        if variables.cleanSOSPoly
            prog = prog.withSDSOS(clean(LNew{sampleIndex}{1},variables.SOSCleanTol));
            prog = prog.withSDSOS(clean(LNew{sampleIndex}{2},variables.SOSCleanTol));
            prog = prog.withSDSOS(clean(-gamma(1)*(xper.'*xper)-DVBal-LNew{sampleIndex}{1}*regionIneq,variables.SOSCleanTol));
            prog = prog.withSDSOS(clean(-gamma(2)*(xper.'*xper)+taudotDenominatorBal - LNew{sampleIndex}{2}*regionIneq,variables.SOSCleanTol));
        else
            prog = prog.withSDSOS(LNew{sampleIndex}{1});
            prog = prog.withSDSOS(LNew{sampleIndex}{2});
            prog = prog.withSDSOS(-gamma(1)*(xper.'*xper)-DVBal-LNew{sampleIndex}{1}*regionIneq);
            prog = prog.withSDSOS(-gamma(2)*(xper.'*xper)+taudotDenominatorBal - LNew{sampleIndex}{2}*regionIneq);
        end
    else
        if variables.cleanSOSPoly
            prog = prog.withSOS(clean(LNew{sampleIndex}{1},variables.SOSCleanTol));
            prog = prog.withSOS(clean(LNew{sampleIndex}{2},variables.SOSCleanTol));
            prog = prog.withSOS(clean(-gamma(1)*(xper.'*xper)-DVBal-LNew{sampleIndex}{1}*regionIneq,variables.SOSCleanTol));
            prog = prog.withSOS(clean(-gamma(2)*(xper.'*xper)+taudotDenominatorBal - LNew{sampleIndex}{2}*regionIneq,variables.SOSCleanTol));
            
        else
            prog = prog.withSOS(LNew{sampleIndex}{1});
            prog = prog.withSOS(LNew{sampleIndex}{2});
            prog = prog.withSOS(-gamma(1)*(xper.'*xper)-DVBal-LNew{sampleIndex}{1}*regionIneq);
            prog = prog.withSOS(-gamma(2)*(xper.'*xper)+taudotDenominatorBal - LNew{sampleIndex}{2}*regionIneq);
        end
    end
    
    opt = spot_sdp_default_options();
    opt.verbose = 0;
    
    if variables.maxGamma_multiplier_step
        solSearchMultiplier = prog.minimize(-sum(gamma),variables.solver,opt);
    else
        solSearchMultiplier = prog.minimize(msspoly(0),variables.solver,opt);
    end
    solAll{sampleIndex} = solSearchMultiplier;
    if ~strcmp(char(solSearchMultiplier.info.status),'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
        display(sprintf('ISSUES FOR SAMPLE %d, %s', sampleIndex, solSearchMultiplier.info.status))
    else
        
    end
    LSolNew{sampleIndex}{1} = solSearchMultiplier.eval(LNew{sampleIndex}{1});
    LSolNew{sampleIndex}{2} = solSearchMultiplier.eval(LNew{sampleIndex}{2});
    
    if variables.findMultiplierBalance
        LSolNew{sampleIndex}{1} = subs(LSolNew{sampleIndex}{1},xper,TiAll{sampleIndex}*xper);
        LSolNew{sampleIndex}{2} = subs(LSolNew{sampleIndex}{2},xper,TiAll{sampleIndex}*xper);
    end
end



end

function [rhoSol, rhodotSol, rhointegralSol, SOSConds] = optimizeRho(precomp, variables, LSol, fTaylorAll, dts, N, TAll, TiAll, opt_rho_Mode, last_Optimal_integral)
xper=variables.xper;
x=variables.x;
uCon=variables.uCon;
tau=variables.tau;

SOSConds = cell(length(precomp),1);

prog = spotsosprog;
prog = prog.withIndeterminate(xper);
prog = prog.withIndeterminate(tau);

% add slack variables if elected
switch variables.rhoStepObjective
    case 'max_single_slack'
        [prog, slack] = prog.newPos(1);
    case 'max_many_slack'
        [prog, slacks] = prog.newPos(40);
end

switch variables.rhoType
    case 'poly_tau'
        tauMono = variables.tauMono;
        % define rho to search over
        [prog, rhoAll] = newFree(prog, variables.rhoDeg + 1);
        rhoPoly = rhoAll'*tauMono;
        rhointegralpoly = integral(rhoPoly,tau);
        rhointegral = subs(rhointegralpoly, tau, precomp{end}.timeNow)  -    subs(rhointegralpoly, tau, 0) ;
        
        iterLastSlide = 40;
    case 'piecewise_linear'
        rhointegral = 0;
        [prog, rhoAll] = newPos(prog,40);
        rhoAll=[rhoAll];
        rhointegral=0;
        iterLastSlide = 39; % last slide is the switching condition...
    case 'chebyshev_poly'
        [prog, rhoAll] = newFree(prog, variables.rhoDeg + 1);
        rhointegral = 0;
        for n = 0:variables.rhoDeg
            tMax = precomp{end}.timeNow;
            tMin = precomp{1}.timeNow;
            
            if n==1
                rhointegrandUpper = (tMax^2)/2;
                rhointegrandLower = (tMin^2)/2;
            else
                rhointegrandUpper = ((n*chebyshevT(n+1,tMax))/(n^2-1) ) - (tMax*chebyshevT(n,tMax))/(n-1) ;
                rhointegrandLower = ((n*chebyshevT(n+1,tMin))/(n^2-1) ) - (tMin*chebyshevT(n,tMin))/(n-1) ;
            end
            rhointegral = rhointegral + rhoAll(n+1)*(rhointegrandUpper - rhointegrandLower);
            iterLastSlide = 40;
        end
end

% set feasibility problem condition: a lower bound for rho integral
switch opt_rho_Mode
    case 'max_rho'
        display('Maximising rho which satisfies transverse Lyapunov criteria...')
        
    case 'rho_feasibility'
        display(sprintf('Recovering strictly feasible rho with last optimised objective (rho_integral) of %f and damping factor of %f ...', last_Optimal_integral, variables.rho_feasibility_step_damping))
        rho_integral_lower_bound = variables.rho_feasibility_step_damping * last_Optimal_integral;
        prog = prog.withPos(rhointegral - rho_integral_lower_bound);
end


for sampleIndex = 1:iterLastSlide
    timeNow = precomp{sampleIndex}.timeNow;
    
    %         SOSConds{sampleIndex} = struct();
    
    switch variables.rhoType
        case 'poly_tau'
            rhoNow = subs(rhoPoly, tau, timeNow);
            rhodot = [];
            
        case 'piecewise_linear'
            rhoNow = rhoAll(sampleIndex);
            rhodot(sampleIndex,1) = (rhoAll(sampleIndex+1)-rhoAll(sampleIndex))/dts(sampleIndex);
            rhointegral = rhointegral+ rhoAll(sampleIndex)*dts(sampleIndex)+.5*rhodot(sampleIndex)*dts(sampleIndex)^2;
        case 'chebyshev_poly'
            rhoNow = chebyshevT(0:variables.rhoDeg, timeNow) * rhoAll;
            rhodot = [];
    end
    
    %% compute transverse conditions...
    [DV, taudotDenominator]=transverseConditions(sampleIndex,precomp{sampleIndex},rhoAll,rhodot,fTaylorAll, variables);
    
    %% enforce slack variables conditions if selected
    switch variables.rhoStepObjective
        case 'max_single_slack'
            prog = prog.withPos(rhoNow - slack);
        case 'max_many_slack'
            prog = prog.withPos(rhoNow - slacks(sampleIndex));
    end
    
    
    DV4All{sampleIndex} = DV;
    DVFull4All{sampleIndex}=DV*taudotDenominator;
    
    % balancing
    if variables.optimizeRhoBalance % Would this actually make sense? how do you recognicile Tdot when you solve everything together...
        % don't bother with coming up with T here as it'd depend on rho
        xperBal = TAll{sampleIndex}*xper;
    else
        xperBal = xper;
    end
    
    DVBal = subs(DV,xper,xperBal);
    L1Bal = subs(LSol{sampleIndex}{1},xper, xperBal);
    L2Bal = subs(LSol{sampleIndex}{2},xper, xperBal);
    taudotDenominatorBal = subs(taudotDenominator, xper, xperBal);
    
    regionIneq = rhoNow - xperBal'*(precomp{sampleIndex}.PMatrixNow)*xperBal;
    
    regionIneq4All{sampleIndex} = regionIneq;
    
    buffer = variables.negativity_buffer_epsilon_rho_step;
    
    if variables.useSDSOS_rho_step
        if variables.cleanSOSPoly
            prog = prog.withSDSOS(clean(-buffer*(xper.'*xper)-DVBal-L1Bal*regionIneq,variables.SOSCleanTol));
            prog = prog.withSDSOS(clean(-buffer*(xper.'*xper)+taudotDenominatorBal - L2Bal*regionIneq,variables.SOSCleanTol));
        else
            prog = prog.withSDSOS(-buffer*(xper.'*xper)-DVBal-L1Bal*regionIneq);
            prog = prog.withSDSOS(-buffer*(xper.'*xper)+taudotDenominatorBal - L2Bal*regionIneq);
        end
    else
        if variables.cleanSOSPoly
            prog = prog.withSOS(clean(-buffer*(xper.'*xper)-DVBal-L1Bal*regionIneq, variables.SOSCleanTol));
            prog = prog.withSOS(clean(-buffer*(xper.'*xper)+taudotDenominatorBal - L2Bal*regionIneq, variables.SOSCleanTol));
        else
            prog = prog.withSOS(-buffer*(xper.'*xper)-DVBal-L1Bal*regionIneq);
            prog = prog.withSOS(-buffer*(xper.'*xper)+taudotDenominatorBal - L2Bal*regionIneq);
        end
    end
    
    % Set some parameters and solve with Sedumi
    clear pars
    % pars.fid = 1;
    pars.cg.maxiter = 20000;
    % pars.cg.refine = 4;
    
    opt = spot_sdp_default_options();
    opt.verbose = 0;
    
end

%% enforce switch condition
sampleIndex = 100;
deltaqDotColumnSpot = variables.deltaqDotColumnSpot;
deltaqDot = reshape(deltaqDotColumnSpot,2,2);
deltaq = variables.deltaq;

switch variables.rhoType
    case 'poly_tau'
        rhoMinus = subs(rhoPoly, tau, precomp{end}.timeNow);
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        rhoPlus = subs(rhoPoly, tau, precomp{1}.timeNow);
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
        
    case 'piecewise_linear'
        
        rhoMinus = rhoAll(end);
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        rhoPlus = rhoAll(1);
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
    case 'chebyshev_poly'
        
        rhoMinus = chebyshevT(0:variables.rhoDeg, precomp{end}.timeNow) * rhoAll;
        V0Minus = precomp{end}.VNow;
        
        xMinus = (precomp{end}.PiNow)'*xper + precomp{end}.xStarNow;
        
        xPlus = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;
        
        xPlusperp = precomp{1}.PiNow*(xPlus - precomp{1}.xStarNow);
        
        rhoPlus = chebyshevT(0:variables.rhoDeg, precomp{1}.timeNow) * rhoAll;
        V0Plus = subs(precomp{1}.VNow,xper,xPlusperp);
        
        
end
% balance
switchingCondition = V0Minus*rhoPlus - V0Plus*rhoMinus ; % needs to be > 0
% balancing
if variables.optimizeRhoBalance % Would this actually make sense? how do you recognicile Tdot when you solve everything together...
    % don't bother with coming up with T here as it'd depend on rho
    xperBal = TAll{sampleIndex}*xper;
else
    xperBal = xper;
end
L1Bal = subs(LSol{sampleIndex}{1},xper, xperBal);
regionIneq = rhoMinus - xperBal'*(precomp{end}.PMatrixNow)*xperBal;
switchingConditionBal = subs(switchingCondition, xper, xperBal);
buffer = variables.negativity_buffer_DISCRETE_rho_step;
% prog = prog.withSOS(V0Minus*rhoPlus - V0Plus*rhoMinus);
if variables.cleanSOSPoly
    if variables.useSDSOS_rho_step
        prog = prog.withSDSOS(clean(-buffer*(xper.'*xper) +switchingConditionBal  -  L1Bal*regionIneq,variables.SOSCleanTol));
    else
        prog = prog.withSOS(clean(-buffer*(xper.'*xper) +switchingConditionBal  -  L1Bal*regionIneq,variables.SOSCleanTol));
    end
else % 'dirty'
    if variables.useSDSOS_rho_step
        prog = prog.withSDSOS(-buffer*(xper.'*xper) +switchingConditionBal  -  L1Bal*regionIneq);
    else
        prog = prog.withSOS(-buffer*(xper.'*xper) +switchingConditionBal  -  L1Bal*regionIneq);
    end
end

if variables.periodicRho_rho_step
    prog = prog.withEqs(rhoPlus - rhoMinus);
end

%% SOS verification...
try
    
    tstart = tic;
    
    switch variables.rhoStepObjective
        case 'max_integral'
            
            %
            switch opt_rho_Mode
                case 'max_rho'
                    solFixMultiplier = prog.minimize(-rhointegral,variables.solver,opt);
                case 'rho_feasibility'
                    %                     solFixMultiplier = prog.minimize(0,           @spot_mosek,opt);
                    solFixMultiplier = prog.minimize(msspoly(0),variables.solver,opt);
            end
        case 'max_single_slack'
            solFixMultiplier = prog.minimize(-slack,variables.solver,opt);
        case 'max_many_slack'
            solFixMultiplier = prog.minimize(-sum(slacks),variables.solver,opt);
    end
    
    solAll = solFixMultiplier;
    if ~strcmp(char(solFixMultiplier.info.status),'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
        display(sprintf('FAILED RHO STEP %s', solFixMultiplier.info.status))
    end
catch err
    
end
SOSConds{end}.sol=solFixMultiplier;
rhoSol = double(solFixMultiplier.eval(rhoAll));
rhodotSol = double(solFixMultiplier.eval(rhodot));
rhointegralSol = double(solFixMultiplier.eval(rhointegral));
program_finished_solve = toc(tstart);

% plot rho
if variables.plot_rho_step
    timeSamples = zeros(40,1);
    rhoResultSamples = zeros(40,1);
    
    for i = 1:40
        timeNow = precomp{i}.timeNow;
        timeSamples(i) = timeNow;
        switch variables.rhoType
            case 'poly_tau'
                rhoPolyResult = rhoSol'*tauMono;
                % gather all the rho...
                rhoResultSamples(i) = subs(rhoPolyResult, tau, timeNow);
            case 'chebyshev_poly'
                rhoResultSamples(i) = chebyshevT(0:variables.rhoDeg, timeNow) * rhoSol;
        end
    end
    rhoSamples = 200;
    rhoSampleFine = zeros(rhoSamples,1);
    sampleTimeVec_rho = linspace(0,precomp{end}.timeNow,rhoSamples);
    for i = 1:rhoSamples
        timeNow = sampleTimeVec_rho(i);
        switch variables.rhoType
            case 'poly_tau'
                rhoPolyResult = rhoSol'*tauMono;
                % gather all the rho...
                rhoSampleFine(i) = subs(rhoPolyResult, tau, timeNow);
            case 'chebyshev_poly'
                rhoSampleFine(i) = chebyshevT(0:variables.rhoDeg, timeNow) * rhoSol;
        end
    end
    figure(); hold on;
    sampled_handle = plot(timeSamples,rhoResultSamples,'x');
    poly_handle = plot(sampleTimeVec_rho, rhoSampleFine);
    title('rho vs sample no.'); xlabel('sample no.'); grid on;
    legend([sampled_handle, poly_handle],'Sampled pts','polynomial')
end
end

function [DV, taudotDenominator] = transverseConditions(sampleIndex,precomp,rhoAll,rhodot,fTaylorAll,variables)
% decompresss here
x = variables.x;
xper = variables.xper;
uCon = variables.uCon;
tau = variables.tau;

timeNow = precomp.timeNow;

dfduNowpp = variables.dfduNowpp;
limCycDuration = variables.limCycDuration;

% get rho and rhodot
switch variables.rhoType
    case 'poly_tau'
        tauMono = variables.tauMono;
        rhoPoly = rhoAll'*tauMono;
        drhoPolydt = diff(rhoPoly,tau);
        rhoNow = subs(rhoPoly,tau,timeNow);
        rhoDotNow = subs(drhoPolydt, tau, timeNow);
    case 'piecewise_linear'
        rhoNow = rhoAll(sampleIndex);
        rhoDotNow = rhodot(sampleIndex);
        
    case 'chebyshev_poly'
        
        rhoNow = chebyshevT(0:variables.rhoDeg, timeNow) * rhoAll;
        
        rhoDotNow = 0;
        for i = 1:variables.rhoDeg
            rhoDotNow = rhoDotNow + rhoAll(i+1) * i*chebyshevU(i-1,timeNow);
        end
        
end
fTaylor = fTaylorAll{sampleIndex};
fNow = precomp.fNow;
PiNow = precomp.PiNow;
PiDotNow = precomp.PiDotNow;
zNow = precomp.zNow;
zDotNow = precomp.zDotNow;
PMatrixNow = precomp.PMatrixNow;
VNow = precomp.VNow;
dVdtauNow = precomp.dVdtauNow;
dVdxperNow = precomp.dVdxperNow;
xStarNow = precomp.xStarNow;
uStarNow = precomp.uStarNow;
RInverse = precomp.RInverse;

% work out u and sub into f
% Get Bperp (transverse dfdu)
dfduNow = dfduNowpp(xStarNow(1),xStarNow(2),xStarNow(3),xStarNow(4));
dtaudu = zNow'*dfduNow/(zNow'*fNow);
transverseB = PiNow*dfduNow - PiNow * fNow * dtaudu;
u = uStarNow + (-1)*RInverse*transverseB'*PMatrixNow*xper;
fTaylorWithFeedback = subs(fTaylor, uCon, u);
fTaylorWithFeedbackxperp = subs(fTaylorWithFeedback, x, xStarNow+PiNow.'*xper);


%% continuous dynamics...
taudotNumerator = zNow.'*fTaylorWithFeedbackxperp;
taudotDenominator = zNow.'*fNow - zDotNow.'*PiNow.'*xper;

% incorporate rho
% recall rhoDot is dRho/dtau
dVdtauNow = dVdtauNow * rhoNow               - VNow * rhoDotNow;   % intentionally skip 1/rho^2 term (to be multiplied in everywhere else)
dVdxperNow = rhoNow  *  dVdxperNow; % 1/rho -> rho because it's multiplied by the rho^2
DV = dVdtauNow*taudotNumerator + dVdxperNow*(taudotNumerator*PiDotNow*PiNow.'*xper ...
    + taudotDenominator*PiNow*fTaylorWithFeedbackxperp - PiNow*fNow*taudotNumerator);
end