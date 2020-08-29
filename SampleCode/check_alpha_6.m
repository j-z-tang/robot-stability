% check taylor series expansion accuracy, check DV, check contraction
clear all

sampleOnSurface = false; % if false, sample anywhere within funnel

% number of samples per surface
samplingTimesPerSurface = 100;


robotParameters = load('alpha_4_lyap_abs_coord_unactuated/z_Pi_P_bezierTheta_extra_stable_ep_1__Q_[100 1 10 0.1]_biggerImpactQ__ZnormalMin_0.06_pOpt100_Qd_factor_1.mat','p');
p = robotParameters.p;

taylorDynamics =  load('alpha_5_taylor_dynamics_abs_unactuated/spline_full_dynamics_taylor_stableEp_1_Q_[100 1 10 0.1]_biggerImpactQ__ZnormalMin_0.06_pOpt_100_Qd_factor_1.mat','p','BtransAll','xperDotTaylorAll','xperDotOriginalDiff','xperDotOriginal','PMatrixAll','QFullNow','PIter','symImpactMap');
xperDotOriginalDiff = taylorDynamics.xperDotOriginalDiff;
xperDotTaylorAll = taylorDynamics.xperDotTaylorAll;
xperDotOriginal = taylorDynamics.xperDotOriginal;
PMatrixAll = taylorDynamics.PMatrixAll;
BtransAll = taylorDynamics.BtransAll;
QFull = taylorDynamics.QFullNow;
symImpactMap  = taylorDynamics.symImpactMap;

verificationLoad = load('alpha_6_ABS_COORD_verification_stableEp_1_Q_[100 1 10 0.1]_ZnormalMin_0.06_pOpt_100_Qd_factor_1/result_5_2_new_z_bilinear_spline_stableEp_1_Q_[100 1 10 0.1]_ZnormalMin_0.06_pOpt_100_SDSOS_1.mat');
fTaylorAll=verificationLoad.fTaylorAll;
f=verificationLoad.f;
dfduNowpp = verificationLoad.dfduNowpp;

% load spot variables first
variables = verificationLoad.variables;
x = variables.x;
xper = variables.xper;
uCon = variables.uCon;
tau = variables.tau;

rhoSolAll = verificationLoad.rhoSolStruct{end};

switch variables.rhoType
    case 'poly_tau'
        rhoPolyResult = rhoSolAll'*variables.tauMono;
        DrhoDtauResult = diff(rhoPolyResult,tau);
    case 'chebyshev_poly'
end

%% sample for DV, taylor series accuracy
deltaqDotColumnSpot = variables.deltaqDotColumnSpot;
deltaqDot = reshape(deltaqDotColumnSpot,2,2);
deltaq = variables.deltaq;

switch variables.rhoType
    case 'poly_tau'
        tauMono = variables.tauMono;
        rhoPlus = subs(rhoSolAll'*tauMono, tau, verificationLoad.precomp{1}.timeNow);
        rhoMinus = subs(rhoSolAll'*tauMono, tau, verificationLoad.precomp{end}.timeNow);
        
    case 'piecewise_linear'
        
    case 'chebyshev_poly'
        
        rhoMinus = chebyshevT(0:variables.rhoDeg, verificationLoad.precomp{end}.timeNow) * rhoSolAll;
        rhoPlus = chebyshevT(0:variables.rhoDeg, verificationLoad.precomp{1}.timeNow) * rhoSolAll;
        
end

V0Minus = verificationLoad.precomp{end}.VNow;

xMinus = (verificationLoad.precomp{end}.PiNow)'*xper + verificationLoad.precomp{end}.xStarNow;

xPlusTaylor = subs(  [deltaq*[x(1:2)]; deltaqDot*[x(3:4)]], x, xMinus) ;

PiNowPlus = verificationLoad.precomp{1}.PiNow;
xStarNowPlus = verificationLoad.precomp{1}.xStarNow;
VNowPlus = verificationLoad.precomp{1}.VNow;

xPlusperpTaylor = verificationLoad.precomp{1}.PiNow*(xPlusTaylor - verificationLoad.precomp{1}.xStarNow);

V0PlusTaylor = subs(verificationLoad.precomp{1}.VNow,xper,xPlusperpTaylor);

for surfaceNo = 40:-1:1
    
    display(sprintf('Checking Lyapunov condition with actual dynamics for %d samples on surface no. %d', samplingTimesPerSurface ,surfaceNo))
    
    precomp = verificationLoad.precomp{surfaceNo};
    fStarNow = precomp.fNow;
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
    timeNow = precomp.timeNow;
%     PdotMatrixNow = precomp.PdotMatrixNow;
    PdotMatrixNow = precomp.PdotMatrixNow;
    
    dfperpdx = xperDotOriginalDiff{surfaceNo};
    xperDot = xperDotOriginal{surfaceNo};
    xperDotTaylor = xperDotTaylorAll{surfaceNo};
    
%     PMatrixNow = PMatrixAll{surfaceNo};
    [U S V] = svd(PMatrixNow);
%     [U S V] = svd(T'*P*T);
    M = U*(sqrt(S));

    taylorSeries = fTaylorAll{surfaceNo};
    dfStardxTaylor = double(  subs(  subs(  diff(taylorSeries,x), x,  xStarNow),  uCon,  uStarNow ) );
    
    % rho and rhodot from results...
    switch variables.rhoType
        case 'poly_tau'
            % gather all the rho...
            rhoNow = double( subs(rhoPolyResult, tau, timeNow) );
            rhoDotNow = double( subs(DrhoDtauResult, tau, timeNow) );
        case 'chebyshev_poly'
            rhoNow = chebyshevT(0:variables.rhoDeg, timeNow) * rhoSolAll;
            rhoDotNow = 0;
            for iter = 1:variables.rhoDeg
                rhoDotNow = rhoDotNow + rhoSolAll(iter+1) * iter*chebyshevU(iter-1,timeNow);
            end
    end
    
    fStarNow =  f(xStarNow,uStarNow);
    dfduNow = dfduNowpp(xStarNow(1),xStarNow(2),xStarNow(3),xStarNow(4));
    
    parfor randomSample = 1:samplingTimesPerSurface
        
        thetaValue = rand()*pi;
        psiValue = rand()*2*pi;
%         xperpSample = (rand(3,1)*sampleRadius*2)-sampleRadius;
        if sampleOnSurface
            r = sqrt(rhoNow);
        else
            r = rand()*sqrt(rhoNow);
        end
        m = [r*sin(thetaValue)*cos(psiValue); r*sin(thetaValue)*sin(psiValue); r*cos(thetaValue)];
        xperpSample = M'\m; % DON'T FORGET TRANSPOSE HERE!!!!!!!
        
        % sanity check
        if sampleOnSurface
            if abs(xperpSample'*PMatrixNow*xperpSample - rhoNow) > 1e-9
                warning('!!!!!!!!!!sampled point does not satisfy equation')
                abs(xperpSample'*PMatrixNow*xperpSample - rhoNow)
                surfaceNo
            end
        else
            if xperpSample'*PMatrixNow*xperpSample > (rhoNow + 1e-9)
                warning('!!!!!!!!!!sampled point does not satisfy equation')
                xperpSample
                xperpSample'*PMatrixNow*xperpSample
                surfaceNo
            end
        end

        xFromxper = xStarNow + PiNow.'*xperpSample;

        dtaudu = zNow'*dfduNow/(zNow'*fStarNow);
        transverseB = PiNow*dfduNow - PiNow * fStarNow * dtaudu;
        %        %u = (-1)*RInverse*transverseB'*PMatrixNow*xper;
        u = uStarNow + (-1)*RInverse*transverseB'*PMatrixNow*xperpSample;
        
        %% get f from taylor series expansion

        
        taylor_eval = double(subs(  subs(taylorSeries, x, xFromxper)  ,uCon,u));
        
        %% evaluate actual (dV/dx)*f from symbolic...
%         f_at_x_eval = f(xFromxper(1),xFromxper(2),xFromxper(3),xFromxper(4),u);
        f_at_x_eval = f(xFromxper,u);
        
        %% compare the f
        if norm(taylor_eval-f_at_x_eval) > 5e-2*norm(f_at_x_eval)
            display('!!!!!!!!!!taylor series f differs from actual f more than tolerance')
            norm(taylor_eval-f_at_x_eval)
        end
        
        %% contraction
        dfdxNow = double(subs(dfperpdx, sym('xper',[3,1]), xperpSample));
        
%         %% compute DV with taylor
        taudotNumerator = zNow.'*taylor_eval;
        taudotDenominator = zNow.'*fStarNow - zDotNow.'*PiNow.'*xperpSample;
        
        dVdtauNowWithRho = dVdtauNow * rhoNow               - VNow * rhoDotNow;   % intentionally skip 1/rho^2 term (to be multiplied in everywhere else)
        dVdxperNowWithRho = rhoNow  *  dVdxperNow; % 1/rho -> rho because it's multiplied by the rho^2
        
        DV_with_taylor =  double(    subs(    dVdtauNowWithRho*taudotNumerator + dVdxperNowWithRho*(taudotNumerator*PiDotNow*PiNow.'*xperpSample ...
            + taudotDenominator*PiNow*taylor_eval - PiNow*fStarNow*taudotNumerator)     ,xper,xperpSample));

        
        
        %% cmopute DV with actual f
        taudotNumerator = zNow.'*f_at_x_eval;
        taudotDenominator = zNow.'*fStarNow - zDotNow.'*PiNow.'*xperpSample;

        dVdtauNowWithRho = dVdtauNow * rhoNow               - VNow * rhoDotNow;   % intentionally skip 1/rho^2 term (to be multiplied in everywhere else)
        dVdxperNowWithRho = rhoNow  *  dVdxperNow; % 1/rho -> rho because it's multiplied by the rho^2

        
        DV_with_actual_f = double(  subs(    dVdtauNowWithRho*taudotNumerator + dVdxperNowWithRho*(taudotNumerator*PiDotNow*PiNow.'*xperpSample ...
            + taudotDenominator*PiNow*f_at_x_eval - PiNow*fStarNow*taudotNumerator)     ,xper,xperpSample)   );
        if norm(DV_with_actual_f - DV_with_taylor)/norm(DV_with_actual_f) > 1e-2
            display(sprintf('!!!!!!!!!DV actal DIFFERENT to DV taylor by > 1 percent on surface No. %d', surfaceNo))
            norm(DV_with_actual_f - DV_with_taylor)/norm(DV_with_actual_f)
        end
        
        if DV_with_actual_f > 0
            display(sprintf('!!!!!!!!!!!DV actual > 0! for surface No.%d', surfaceNo))
            DV_with_actual_f
        end
        
        if DV_with_taylor >0
            display(sprintf('!!!!!!!!!!DV taylor > 0! for surfaceNo. %d',surfaceNo))
            DV_with_taylor
        end
        
        dtaudx = (zNow'*dfStardxTaylor*PiNow' + zDotNow'*PiNow')/(zNow'*fStarNow);
        transLinA = (PiDotNow*PiNow' + PiNow*dfStardxTaylor*PiNow' - PiNow*fStarNow*dtaudx);
        
       %% checking switching condition
       if surfaceNo == 40
           if randomSample == 1; display('Checking switching condition'); end
           

            switchingConditionTaylor = double(  subs( V0Minus*rhoPlus - V0PlusTaylor*rhoMinus  ,xper, xperpSample)) ; % needs to be > 0 
           
            if switchingConditionTaylor < -1e-10  % floating point accuracy
                display('!!!!!!!!!!!!SWITCHING CONDITION SAMPLE VIOLATED (Taylor expansion of delta)')
                display(xperpSample)
                display(switchingConditionTaylor)
            end
            
            [deltaqActual, deltaqDotActual] = symImpactMap(xFromxper,p);
            xPlus = subs([deltaqActual*[x(1:2)]; deltaqDotActual*[x(3:4)]], x, xMinus) ;
            
            % check how accurate is taylor series expansion of impact map
            deltaqDotTaylor = double( subs(deltaqDot, x, xFromxper)  );
            deltaqDot_difference = deltaqDotActual - deltaqDotTaylor;
            
            if norm(deltaqDot_difference)/norm(deltaqDotActual) > 1e-2
                
                display('!!!! taylor series of impact map differs from the actual impact map by more than 1 percent')
                
            end
            
            
            xPlusperp = PiNowPlus*(xPlus - xStarNowPlus);
            V0Plus_actual = subs(VNowPlus, xper, xPlusperp);
            
            switchingConditionActual = double(  subs( V0Minus*rhoPlus - V0Plus_actual*rhoMinus  ,xper, xperpSample)) ; % needs to be > 0 
            if switchingConditionActual < -1e-10  % floating point accuracy
                display('!!!!!!!!!!!!SWITCHING CONDITION SAMPLE VIOLATED (actual)')
                display(xperpSample)
                display(switchingConditionActual)
            end
            
       end
        
    end
end
