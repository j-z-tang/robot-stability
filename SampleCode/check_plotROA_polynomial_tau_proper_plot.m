% plot the ROA with proper projections and shading etc.

QFullListSize = 1;
actuated = false;
rel_coord = false;
PFL = false;

% loadFolder = 'alpha_5_taylor_dynamics'
% loadFolder = 'Unactuated_alpha_5_taylor_dynamics'

if actuated
    %     savePlotFolder = 'ROA_plots'
    %
    if rel_coord
        if PFL
            savePlotFolder = 'ROA_plots_PFL'
        else
            savePlotFolder = 'ROA_plots_rel_coord'
        end
    else %abs coord
        if PFL
            savePlotFolder = 'ROA_plots_abs_coord_PFL'
        else
            savePlotFolder = 'ROA_plot_abs_coord_full'
        end
    end
    
else
    savePlotFolder = 'ROA_plots_abs_coord_nominal_traj_energy_optimal';
end
if ~exist(savePlotFolder,'dir')
    mkdir(savePlotFolder)
end

plot_random_samples_of_xperp = false;
plot_z_dir = false;

plot_nominal_traj = false;

funnel_whiteness = 0.9;
funnel_shell_visibility = 0.05;
fainter_outline = true;
edge_color = 0.5;

set_axis_ROA_plot = true;
axis_ROA_plot = [-0.5,0.5,-2.5,3.5]; % for abs coordinates

for SDSOS = [1]
for normalDeltaNow=[0.06]
% for stableEpsilonNow = [0.1]
for stableEpsilonNow = [1]
% for discreteQFactor = [1]
for discreteQFactor = [1]
% for QFullIndex=1:length(QFullList)
for QFullIndex=1:QFullListSize
% QFullList={diag([100 1 10 0.1]), diag([50 1 5 0.5]), diag([1 1 1 1]) ,diag([1 100 0.1 10]) };
QFullList={diag([100 1 10 0.1])};
% QFullList={diag([1 1 1 1])};
% QFullList={diag([100 1 10 0.1]), diag([50 1 5 0.5]) }
QFullNow=QFullList{QFullIndex};
% if QFullIndex == 2
%     stableEpsilonNow = 0;
% end
pOpt = 100;

useVerificationSol = true; % are we using results from verification? Set 0 if just seeing shape of P

plotRhoDot = false;
plotDirectionArrowsForPhase = false;

fileName = sprintf('result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_SDSOS_%d.mat',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,SDSOS);
if actuated
%     if SDSOS
%         loadFolderMaster = 'alpha_6_verification_SDSOS_ALL_MOSEK';
%     else
%         loadFolderMaster = 'alpha_6_verification_SOS_ALL_MOSEK';
%     end
%     loadFolderName = sprintf('alpha_6_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
%     fileName = sprintf('result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d.mat',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt);
    %     filePath = sprintf('%s/%s/%s',loadFolderMaster, loadFolderName, fileName)
    if rel_coord
        if PFL
            loadFolderName = sprintf('alpha_6_PFL_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
        else % rel coord full dynamics
            loadFolderName = sprintf('alpha_6_REL_COORD_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
        end


    else % abs_coord
        if PFL
            loadFolderName = sprintf('alpha_6_PFL_abs_coord_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
        else 
            loadFolderName = sprintf('alpha_6_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
            fileName = sprintf('result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d.mat',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt);
%             loadFolderName = sprintf('alpha_6_ABS_COORD_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
%             fileName = sprintf('result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_SDSOS_%d.mat',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,SDSOS);
        end
    end
    
    filePath = sprintf('%s/%s',loadFolderName, fileName);
    
else
%     loadFolderName = sprintf('Unactuated_alpha_6_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);
%     fileName = sprintf('result_5_2_new_z_bilinear_spline_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_SDSOS_%d.mat',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,SDSOS);
%     filePath = sprintf('%s/%s',loadFolderName, fileName)
        loadFolderName = sprintf('alpha_6_ABS_COORD_verification_stableEp_%d_Q_%s_ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d',stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor);  

        filePath = sprintf('%s/%s',loadFolderName, fileName);
end



if ~exist(filePath,'file')
    display('NO RESULTS FOR THIS GUY')
    continue
end

if useVerificationSol
%     load result_5_2_new_z_bilinear_spline.mat precomp rhoSolStruct Pipp zpp f zDotpp TAllStruct
%     load invariance/result_5_2_bilinear_PARTIAL.mat precomp rhoSolStruct Pipp zpp f zDotpp TAllStruct
%        load result_5_2_bilinear_PARTIAL precomp rhoSolStruct Pipp zpp f zDotpp TAllStruct
%       load alpha_6_verification/result_5_2_new_z_bilinear_spline_stableEp_0_Q_[50 1 5 0.5]_ZnormalMin_0.01_pOpt_100.mat precomp variables rhoSolStruct Pipp zpp f zDotpp TAllStruct
%         load alpha_6_verification/result_5_2_new_z_bilinear_spline_stableEp_1_Q_[100 1 10 0.1]_ZnormalMin_0.01_pOpt_100.mat ... 
% load Unactuated_alpha_6_verification/result_5_2_new_z_bilinear_spline_stableEp_1_Q_[100 1 10 0.1]_ZnormalMin_0.06_pOpt_100.mat ... 
%             precomp variables xStarpp nominal_u rhoSolStruct Pipp zpp f zDotpp TAllStruct stableEpsilonNow QFull sampleStance sampleTimeVec
load(filePath,'multiplier_step_rho_damping_all','rho_feasibility_damping_all','rhointegralAll', 'precomp', 'variables', 'xStarpp', 'nominal_u', 'rhoSolStruct', 'Pipp', 'zpp', 'f','fTaylorAll', 'zDotpp', 'TAllStruct', 'stableEpsilonNow', 'QFull', 'sampleStance', 'sampleTimeVec','dfduNowpp');
        
else 
%     load z_Pi_P_splines_extra_stable_ep_1__Q_10_5_1_0.5.mat PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull stateAll limCycDuration
%     load z_Pi_P_splines_extra_stable_ep_2__Q_5_1_1_0.5.mat PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull stateAll limCycDuration
%     load z_Pi_P_splines_extra_stable_ep_1__Q_100_50_10_5__ZnormalMin_0.10_ZsmoothCost_10_dz2dtCoef_2_dotProdCoef_96.mat PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull stateAll limCycDuration
%     load pre_2_1_mat_hipMass50/z_Pi_P_splines_extra_stable_ep_0__Q_[100 1 10 0.1]_biggerImpactQ__ZnormalMin_0.10_ZsmoothCost_85_dz2dtCoef_48_dotProdCoef_8.mat ...
% load pre_2_1_mat/z_Pi_P_bezierTheta_extra_stable_ep_0__Q_[1000 1 10 0.01]_biggerImpactQ__ZnormalMin_0.00_pOpt200.mat ...
%     load pre_2_1_mat/z_Pi_P_bezierTheta_extra_stable_ep_1__Q_[100 1 10 0.1]_biggerImpactQ__ZnormalMin_0.00_pOpt200.mat ...
%         PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull stateAll limCycDuration stableEpsilon smoothCostCoef dz2dtCoef dotProdCoef
%     sampleNo = 40;
%     sampleTimeVec = 0:limCycDuration/(sampleNo-1):limCycDuration;
end






subplotX = 4;
subplotY = 2;

sampleNo = 40;
% sampleTimeVec = 0:stancePeriod/(sampleNo-1):stancePeriod;
zDirHandle=[];
colourMap = winter(sampleNo);

swingAxMin = -0.5;
swingAxMax = 0.5;
stanceAxMin = -0.45;
stanceAxMax = 0.4;

arrowSize = 0.1;
switchSurfSize = 0.2;
% maxNormSurfSize = 0.139;
% maxNormSurfSize = 0.1163;
maxNormSurfSize = 1.5;
for i=1:sampleNo
    sampleTime = sampleTimeVec(i);
    xStarNow = xStarpp(sampleTime);
    uStarNow = nominal_u(sampleTime);
    qSwAll(i,1)=xStarNow(1);
    qStAll(i,1)=xStarNow(2);
    qSwDotAll(i,1)=xStarNow(3);
    qStDotAll(i,1)=xStarNow(4);
    
    xStarAll(:,i) = xStarNow;
    uStarAll(:,i) = uStarNow;
    zppAll(:,i) = zpp(sampleTime); % each z a column
    zDotppAll(:,i) = zDotpp(sampleTime);
    fAll(:,i) = f(xStarNow,uStarNow);
end








% DONE: check why the middle of the cycle 'radius' of the ellipse, aka rho, when plotting
radius = .03;

% plot dimension
xdimension = 4; % dimension of the original variable x
plotDimension = {[1,3],[2,4]}; % the dimensions to be plotted

tau = variables.tau;
x = variables.x;
xper = variables.xper;
uCon = variables.uCon;

%% plot limit cycle
sampleNo = 40;
% fetch the limit cycle from precomp...
for i=1:sampleNo
    if useVerificationSol
        sampleTime = precomp{i}.timeNow;
        xStarNow = precomp{i}.xStarNow;
        qSwAll(i,1)=xStarNow(1);
        qStAll(i,1)=xStarNow(2);
        qSwDotAll(i,1)=xStarNow(3);
        qStDotAll(i,1)=xStarNow(4);
    else
        sampleTime = sampleTimeVec(i);
        xStarNow = xStarpp(sampleTime);
        qSwAll(i,1)=xStarNow(1);
        qStAll(i,1)=xStarNow(2);
        qSwDotAll(i,1)=xStarNow(3);
        qStDotAll(i,1)=xStarNow(4);
        
    end
    
    zppAll(:,i) = zpp(sampleTime); % each z a column
    zDotppAll(:,i) = zDotpp(sampleTime);
%     fAll(:,i) = f(xStarNow);
end
% sample more finely for xStar to plot...
fine_sample_times = 200;
final_time = precomp{end}.timeNow;
sampleTimeFineAll = linspace(0,final_time,fine_sample_times);

for i = 1:fine_sample_times
    finetimeNow = sampleTimeFineAll(i);
    xStarNow = xStarpp(finetimeNow);
    qSwAllFine(i,1)=xStarNow(1);
    qStAllFine(i,1)=xStarNow(2);
    qSwDotAllFine(i,1)=xStarNow(3);
    qStDotAllFine(i,1)=xStarNow(4);
end

%% plot transversal surfaces
h=[];
colors = lines(40);
samplexperTimes_ForFunnel = 1000;
samplexperTimes_ForVDot = 1000;
if useVerificationSol
    rhoSolAll = rhoSolStruct{end};
    TMatrices=TAllStruct{end};
else
%     rhoSolAll = radius*ones(40,1);
end

switch variables.rhoType
    case 'poly_tau'
        rhoPolyResult = rhoSolAll'*variables.tauMono;
        DrhoDtauResult = diff(rhoPolyResult,tau);
    case 'chebyshev_poly'
end















%% plot DV
% checkDVNegative = false;
% sampleOnSurface = false;
% minSampleRhoFactor = 1.0;
% maxSampleRhoFactor = 2.0;
% 
% minDVAllSurf = zeros(1,40);
% maxDVAllSurf = zeros(1,40);
% verifiedRho = zeros(1,40);
% 
% for surfaceNo = 1:40
%     surfaceNo
%     fStarNow = precomp{surfaceNo}.fNow;
%     PiNow = precomp{surfaceNo}.PiNow;
%     PiDotNow = precomp{surfaceNo}.PiDotNow;
%     zNow = precomp{surfaceNo}.zNow;
%     zDotNow = precomp{surfaceNo}.zDotNow;
%     PMatrixNow = precomp{surfaceNo}.PMatrixNow;
%     VNow = precomp{surfaceNo}.VNow;
%     dVdtauNow = precomp{surfaceNo}.dVdtauNow;
%     dVdxperNow = precomp{surfaceNo}.dVdxperNow;
%     xStarNow = precomp{surfaceNo}.xStarNow;
%     uStarNow = precomp{surfaceNo}.uStarNow;
%     RInverse = precomp{surfaceNo}.RInverse;
%     timeNow = precomp{surfaceNo}.timeNow;
% %     PdotMatrixNow = precomp.PdotMatrixNow;
%     PdotMatrixNow = precomp{surfaceNo}.PdotMatrixNow;
%     
%     [U S V] = svd(PMatrixNow);
% %     [U S V] = svd(T'*P*T);
%     M = U*(sqrt(S));
% 
%     taylorSeries = fTaylorAll{surfaceNo};
%     dfStardxTaylor = double(  subs(  subs(  diff(taylorSeries,x), x,  xStarNow),  uCon,  uStarNow ) );
%     
%      % rho and rhodot from results...
%     switch variables.rhoType
%         case 'poly_tau'
%             % gather all the rho...
%             rhoNow = double( subs(rhoPolyResult, tau, timeNow) );
%             rhoDotNow = double( subs(DrhoDtauResult, tau, timeNow) );
%         case 'chebyshev_poly'
%             rhoNow = chebyshevT(0:variables.rhoDeg, timeNow) * rhoSolAll;
%             rhoDotNow = 0;
%             for iter = 1:variables.rhoDeg
%                 rhoDotNow = rhoDotNow + rhoSolAll(iter+1) * iter*chebyshevU(iter-1,timeNow);
%             end
%     end
%     verifiedRho(surfaceNo) = rhoNow;
%     
%     fStarNow =  f(xStarNow,uStarNow);
%     dfduNow = dfduNowpp(xStarNow(1),xStarNow(2),xStarNow(3),xStarNow(4));
%     
%     xSamples_THIS_SURF = zeros(4,samplexperTimes_ForVDot);
%     DV_with_actual_f_THIS_SURF = zeros(1,samplexperTimes_ForVDot);
%     Positive_DV_rho = zeros(1,samplexperTimes_ForVDot);
%     
%     parfor randomSample = 1:samplexperTimes_ForVDot
%         
%         thetaValue = rand()*pi;
%         psiValue = rand()*2*pi;
% %         xperpSample = (rand(3,1)*sampleRadius*2)-sampleRadius;
%         if sampleOnSurface
%             r = sqrt(rhoNow);
%         else
%             r = ((maxSampleRhoFactor-minSampleRhoFactor)*rand() + minSampleRhoFactor) *sqrt(rhoNow);
%         end
%         m = [r*sin(thetaValue)*cos(psiValue); r*sin(thetaValue)*sin(psiValue); r*cos(thetaValue)];
%         xperpSample = M'\m; % DON'T FORGET TRANSPOSE HERE!!!!!!!
%         
%         % sanity check
%         if sampleOnSurface
%             if abs(xperpSample'*PMatrixNow*xperpSample - rhoNow) > 1e-9
%                 warning('!!!!!!!!!!sampled point does not satisfy equation')
%                 abs(xperpSample'*PMatrixNow*xperpSample - rhoNow)
%                 surfaceNo
%             end
%         else
%             if xperpSample'*PMatrixNow*xperpSample > (r^2 + 1e-9)
%                 warning('!!!!!!!!!!sampled point does not satisfy equation')
%                 xperpSample
%                 xperpSample'*PMatrixNow*xperpSample
%                 surfaceNo
%             end
%         end
% 
%         xFromxper = xStarNow + PiNow.'*xperpSample;
%         xSamples_THIS_SURF(:,randomSample) = xFromxper;
%          
%         
% %         fStarNow =  f(xStarNow(1), xStarNow(2), xStarNow(3), xStarNow(4), uStarNow);
%         %% compute controller here
% 
% %         dfduNow =[   0
% %             0
% %             (4*(2*cos(qSt - qSw) - 13))/(5*(4*cos(qSt - qSw)^2 - 13))
% %             -(4*(2*cos(qSt - qSw) - 1))/(5*(4*cos(qSt - qSw)^2 - 13))];
% %         if PFL
% %             dfduNow = dfdu_PFL_rel(xStarNow,p);
% %         end
% 
%         dtaudu = zNow'*dfduNow/(zNow'*fStarNow);
%         transverseB = PiNow*dfduNow - PiNow * fStarNow * dtaudu;
%         %        %u = (-1)*RInverse*transverseB'*PMatrixNow*xper;
%         u = uStarNow + (-1)*RInverse*transverseB'*PMatrixNow*xperpSample;
%         
%         %% get f from taylor series expansion
% 
%         
% %         taylor_eval = double(subs(  subs(taylorSeries, x, xFromxper)  ,uCon,u));
%         
%         %% evaluate actual (dV/dx)*f from symbolic...
% %         f_at_x_eval = f(xFromxper(1),xFromxper(2),xFromxper(3),xFromxper(4),u);
%         f_at_x_eval = f(xFromxper,u);
%         
% %         %% compare the f
% %         if norm(taylor_eval-f_at_x_eval) > 5e-2*norm(f_at_x_eval)
% %             display('!!!!!!!!!!taylor series f differs from actual f more than tolerance')
% %             norm(taylor_eval-f_at_x_eval)
% %         end
%         
% 
% % %         %% compute DV with taylor
% %         taudotNumerator = zNow.'*taylor_eval;
% %         taudotDenominator = zNow.'*fStarNow - zDotNow.'*PiNow.'*xperpSample;
% %         
% % %         dVdtauNowWithRho = taudotDenominator*(dVdtauNow * rhoNow)  -  taudotNumerator*(VNow * rhoDotNow);
% % %         dVdxperNowWithRho = taudotDenominator*(rhoNow  *  dVdxperNow);
% %         
% %         dVdtauNowWithRho = dVdtauNow * rhoNow               - VNow * rhoDotNow;   % intentionally skip 1/rho^2 term (to be multiplied in everywhere else)
% %         dVdxperNowWithRho = rhoNow  *  dVdxperNow; % 1/rho -> rho because it's multiplied by the rho^2
% %         
% %         DV_with_taylor =  double(    subs(    dVdtauNowWithRho*taudotNumerator + dVdxperNowWithRho*(taudotNumerator*PiDotNow*PiNow.'*xperpSample ...
% %             + taudotDenominator*PiNow*taylor_eval - PiNow*fStarNow*taudotNumerator)     ,xper,xperpSample));
% 
%         %% cmopute DV with actual f
%         taudotNumerator = zNow.'*f_at_x_eval;
%         taudotDenominator = zNow.'*fStarNow - zDotNow.'*PiNow.'*xperpSample;
%         
% %         dVdtauNowWithRho = taudotDenominator*(dVdtauNow * rhoNow)  -  taudotNumerator*(VNow * rhoDotNow);
% %         dVdxperNowWithRho = taudotDenominator*(rhoNow  *  dVdxperNow);
% 
%         dVdtauNowWithRho = dVdtauNow * rhoNow               - VNow * rhoDotNow;   % intentionally skip 1/rho^2 term (to be multiplied in everywhere else)
%         dVdxperNowWithRho = rhoNow  *  dVdxperNow; % 1/rho -> rho because it's multiplied by the rho^2
% 
%         
%         DV_with_actual_f = double(  subs(    dVdtauNowWithRho*taudotNumerator + dVdxperNowWithRho*(taudotNumerator*PiDotNow*PiNow.'*xperpSample ...
%             + taudotDenominator*PiNow*f_at_x_eval - PiNow*fStarNow*taudotNumerator)     ,xper,xperpSample)   );
% %         if norm(DV_with_actual_f - DV_with_taylor)/norm(DV_with_actual_f) > 1e-2
% %             display(sprintf('!!!!!!!!!DV actal DIFFERENT to DV taylor by > 1 percent on surface No. %d', surfaceNo))
% %             norm(DV_with_actual_f - DV_with_taylor)/norm(DV_with_actual_f)
% %         end
%         if checkDVNegative
%             if DV_with_actual_f > 0
%                 display(sprintf('!!!!!!!!!!!DV actual > 0! for surface No.%d', surfaceNo))
%                 DV_with_actual_f
%             end
%         end
%         if DV_with_actual_f > 0
%             Positive_DV_rho(randomSample) = r^2;
%         end
% %         if DV_with_taylor >0
% %             display(sprintf('!!!!!!!!!!DV taylor > 0! for surfaceNo. %d',surfaceNo))
% %             DV_with_taylor
% %         end
%         DV_with_actual_f_THIS_SURF(randomSample) = DV_with_actual_f;
%         
%     end
%     
%     DV_with_actual_f_ALL_SURFACES{surfaceNo} = DV_with_actual_f_THIS_SURF;
%     minDVAllSurf(surfaceNo) = min(DV_with_actual_f_THIS_SURF);
%     maxDVAllSurf(surfaceNo) = max(DV_with_actual_f_THIS_SURF);
%     xSamples_ALL_SURFACESW{surfaceNo} = xSamples_THIS_SURF;
%     Positive_DV_rho_ALL_SURFACES{surfaceNo} = Positive_DV_rho;
% end
% minDVValue = min(minDVAllSurf); %largest negative
% maxDVValue = max(maxDVAllSurf); % smallest negative
% save(sprintf('%s/Sample_DV_data_Q_%s_ep_%d_Qd_factor_%d_SDSOS=%d.mat',savePlotFolder,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,SDSOS)); 
% load(sprintf('%s/Sample_DV_data_Q_%s_ep_%d_Qd_factor_%d_SDSOS=%d.mat',savePlotFolder,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,SDSOS)); 




















% set up figure
rho_ROA_fig_Handle = figure('PaperType','a3'); 
orient(rho_ROA_fig_Handle,'portrait')
paperUnits = get(rho_ROA_fig_Handle, 'PaperUnits');
set(rho_ROA_fig_Handle,'PaperUnits','centimeters');
paperSize = get(rho_ROA_fig_Handle,'PaperSize');
paperPosition = [1.5 1.5 paperSize - 1.5];
set(rho_ROA_fig_Handle,'PaperPosition', paperPosition);
set(rho_ROA_fig_Handle,'PaperUnits',paperUnits);
colors = lines(40);
% set position of graph on screen so we can see what's happenin
widthPixels = 900;
heightPixels = 1200;
f99Pos = get(rho_ROA_fig_Handle,'Position');
f99Pos = [f99Pos(1), f99Pos(2)-240, widthPixels, heightPixels];
set(rho_ROA_fig_Handle, 'Position', f99Pos)

hold on;











%% plot Funnel
% plot the swing leg
figure(rho_ROA_fig_Handle);
hold on
xlabel('$\theta$ (rads)', 'interpreter','LaTeX');
ylabel('$\dot\theta$ (rads/s)', 'interpreter','LaTeX');
% title('\rho = 0.005 constant with P generated with \epsilon = 0, Q=diag([10 5 1 0.5])')
if useVerificationSol
    title(sprintf('ROA from SDSOS w/ ep=%d, Q=%s', stableEpsilonNow, mat2str(diag(QFull)')))
%     title('outline of "invariant funnel" from transverse dynamics and SOS')
else
    title(sprintf('rho=%.3f w/ smooth=%d, dz2dt=%d, dotProd=%d ep=%d, Q=%s',radius, smoothCostCoef, dz2dtCoef, dotProdCoef, stableEpsilon, mat2str(diag(QFull)')))
end
qSwPlotHandle = plot(qSwAll,qSwDotAll);
qStPlotHandle = plot(qStAll,qStDotAll);
if plot_nominal_traj;
qSwPlotHandle = plot(qSwAllFine,qSwDotAllFine,'b:','LineWidth',2);
qStPlotHandle = plot(qStAllFine,qStDotAllFine,'b:','LineWidth',2);
end
drawnow


%% iterate through each sample point
for i = 1:40
    figure(rho_ROA_fig_Handle);
    
    if set_axis_ROA_plot
        axis(axis_ROA_plot);
    end
    
    % for i=1
    if useVerificationSol
        timeNow = precomp{i}.timeNow;
        Pi = precomp{i}.PiNow;
        P = precomp{i}.PMatrixNow;
        %     P=eye(3);
        x0 = precomp{i}.xStarNow;
%         T=TMatrices{i};
    else 
%         t = sampleTimeVec(i);
%         Pi = Pipp(t);
%         x0 = xStarpp(t);
%         P = PMatrixpp(t);
    end
    

    switch variables.rhoType
        case 'poly_tau'
            rhoPolyResult = rhoSolAll'*variables.tauMono;
            % gather all the rho...
            rhoValue = double(subs(rhoPolyResult, tau, timeNow));
        case 'chebyshev_poly'
            rhoValue = chebyshevT(0:variables.rhoDeg, timeNow) * rhoSolAll;
    end
    
    [U S V] = svd(P);
%     [U S V] = svd(T'*P*T);
    M = U*(sqrt(S));
    
    % sample xper and plot
    samplex_this_surf =  zeros(4,samplexperTimes_ForFunnel);
    parfor j=1:samplexperTimes_ForFunnel
        thetaValue = rand()*pi;
        psiValue = rand()*2*pi;
        
         %% use symbolic maths to solve for x that lies on the transverse ellipsoid
%         rSolAll = subs(subs( subs(subs(solr, rho, rhoValue), theta, thetaValue), psi, psiValue), Psym, P);
%         
%         rSol = rSolAll(rSolAll>0); % only include r that are bigger than zero
%         
%         xperSol = double(subs(subs(subs(xper,theta,thetaValue),psi,psiValue),r, rSol)); % sub r, theta, psi in xper
%         
%         samplexNow = x0 + Pi'*xperSol; % in global frame
%         samplex{i}(:,j) = samplexNow; % each sample is a column

        %% use SVD
        % note xper = inv(M)*m; also m*inv(M)*P*inv(M)'*m = r  <--->  m*I*m = r
        r = sqrt(rhoValue);
        m = [r*sin(thetaValue)*cos(psiValue); r*sin(thetaValue)*sin(psiValue); r*cos(thetaValue)];
        xperSol = M'\m; % DON'T FORGET TRANSPOSE HERE!!!!!!!
        
        % sanity check
        if abs(xperSol'*P*xperSol - rhoValue) > 1e-9
            warning('sampled point does not satisfy equation')
            abs(xperSol'*P*xperSol - rhoValue)
            i
        end
        
        samplexNow = x0 + Pi'*xperSol; % in global frame
        samplex_this_surf(:,j) = samplexNow; % each sample is a column
        
        % plot the sample
        if plot_random_samples_of_xperp
            for k = 1:length(plotDimension)
                plotdims = plotDimension{k};
                plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'.','Color',colors(i,:));
            end
        end
    end
    samplex{i} = samplex_this_surf; % each sample is a column
    
    
    
    
    %% Shade the funnel
    
    for k = 1:length(plotDimension)
        plotdims = plotDimension{k};
        
        ellipse_conv_Hull = convhull(  samplex_this_surf(plotdims(1),:), samplex_this_surf(plotdims(2),:)  );
        
        ellipse_outline_X = samplex_this_surf(plotdims(1),ellipse_conv_Hull);
        ellipse_outline_Y = samplex_this_surf(plotdims(2),ellipse_conv_Hull);
        
        
        
        % fill in funnel
        if i~=1
            
            funnel_boundary_samples_X = [ ellipse_outline_X , last_ellipse_outline_X{k} ];
            funnel_boundary_samples_Y = [ ellipse_outline_Y , last_ellipse_outline_Y{k} ];
            
            funnel_boundary_convex_hull = convhull(funnel_boundary_samples_X, funnel_boundary_samples_Y);
            
%             funnel_fill_handle = fill( samplex_this_surf(plotdims(1),ellipse_conv_Hull), samplex_this_surf(plotdims(2),ellipse_conv_Hull) ,   'b');
            funnel_fill_handle = fill( funnel_boundary_samples_X(funnel_boundary_convex_hull),  funnel_boundary_samples_Y(funnel_boundary_convex_hull), 'b');

            set(funnel_fill_handle,'FaceColor',  [funnel_whiteness,funnel_whiteness,funnel_whiteness])
            if fainter_outline; 
                set(funnel_fill_handle,'EdgeColor',[edge_color,edge_color,edge_color])
            end
            set(funnel_fill_handle,'facealpha',funnel_shell_visibility)
        end
        
        ellipse_outline_handle = plot(   ellipse_outline_X,   ellipse_outline_Y   ,'r');
        % fill in last guy
        if i == 40
            funnel_fill_handle = fill(  ellipse_outline_X,   ellipse_outline_Y   ,'r');
            set(funnel_fill_handle,'FaceColor',  [funnel_whiteness,funnel_whiteness,funnel_whiteness])
            if fainter_outline
                set(funnel_fill_handle,'EdgeColor',  [edge_color,edge_color,edge_color])
            end
            set(funnel_fill_handle,'facealpha',funnel_shell_visibility)
        end
        
        last_ellipse_outline_X{k} = ellipse_outline_X;
        last_ellipse_outline_Y{k} = ellipse_outline_Y;
        
    end
    
    
    
    
    
    

    
    
    
    %% plot z
    if plot_z_dir
        zNow = zppAll(:,i);
        z1 = x0;
        z2 = x0+(zNow*0.2);
        for j = 1:length(plotDimension)
            plotdims = plotDimension{j};
            h=[h,line([z1(plotdims(1));z2(plotdims(1))],[z1(plotdims(2));z2(plotdims(2))],'Color',[1 0 0])];
        end
    end
    drawnow
end









%% plot controller trajectory
% % controllerLoaded = load('controller_sim_slow_fast_damp_0.829250.mat');
% % controllerLoaded = load('controller_sim_slow_fast_damp_0.829250_WORKS_COMPARE[100_10_1_0.1]_WITH_WITHOUT_OFFSET.mat');
% % 
% 
% controllerLoaded = load('controller_sim_slow_fast_damp_0.830000.mat');
% slowControllerTraj_XFULL = controllerLoaded.slowControllerTraj_XFULL;
% fastControllerTraj_XFULL = controllerLoaded.fastControllerTraj_XFULL;
% 
% 
% qSwAllSlow = slowControllerTraj_XFULL(:,1);
% qSwDotAllSlow = slowControllerTraj_XFULL(:,3);
% qStAllSlow = slowControllerTraj_XFULL(:,2);
% qStDotAllSlow = slowControllerTraj_XFULL(:,4);
% 
% 
% initial_cond_handle = plot(qSwAllSlow(1), qSwDotAllSlow(1), 'go','LineWidth',2);
% initial_cond_handle = plot(qStAllSlow(1), qStDotAllSlow(1), 'go','LineWidth',2);
% 
% qSwPlotHandle_slow = plot(qSwAllSlow,qSwDotAllSlow,'b-.','LineWidth',1.5);
% qStPlotHandle_slow = plot(qStAllSlow,qStDotAllSlow,'b-.','LineWidth',1.5);
% 
% for steps=1:3
% qSwAllFast = fastControllerTraj_XFULL{steps}(:,1);
% qSwDotAllFast =  fastControllerTraj_XFULL{steps}(:,3);
% qStAllFast = fastControllerTraj_XFULL{steps}(:,2);
% qStDotAllFast =  fastControllerTraj_XFULL{steps}(:,4);
% 
% qSwPlotHandle_fast = plot(qSwAllFast,qSwDotAllFast,'k','LineWidth',1.5);
% qStPlotHandle_fast = plot(qStAllFast,qStDotAllFast,'k','LineWidth',1.5);
% end
% if plot_nominal_traj
%     legend([initial_cond_handle, qSwPlotHandle_slow,qSwPlotHandle_fast, qSwPlotHandle], 'Initial Condition', 'Trajectory for controller with small stability region' , 'Trajectory for controller with large stability region' ,'Nominal virtual constraint trajectory')
% else
%     legend([initial_cond_handle, qSwPlotHandle_slow,qSwPlotHandle_fast], 'Initial Condition', 'Trajectory for controller with small stability region' , 'Trajectory for controller with large stability region')
%     
% end








%% plot the DV
% figure(rho_ROA_fig_Handle)
% for i = 1:40
%     
%     DVThisSurf = DV_with_actual_f_ALL_SURFACES{i};
%     xSamplesThisSurf = xSamples_ALL_SURFACESW{i};
%     xStarNow = precomp{i}.xStarNow;
%     for j = 1:length(DVThisSurf)
%         samplexNow = xSamplesThisSurf(:,j);
%         DVNow = DVThisSurf(j);
%         for k = 1:length(plotDimension)
%             plotdims = plotDimension{k};
%             if DVNow < 0 % DV Negative
% %                 if abs(DVNow) < 0.1*abs(minDVValue)
% %                     plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'.','Color',[1-abs(DVNow/(0.2*minDVValue)),1-abs(DVNow/(0.2*minDVValue)),0]);
% %                 else
% %                     plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'.','Color',[0,0,0]);
% %                 end
%             else % DV is positive
% %                 plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'x','Color',[1,0,0]);
%                 plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'x','Color',colors(i,:));
%                 plot(xStarNow(plotdims(1)),xStarNow(plotdims(2)),'*','Color',colors(i,:));
%                 plot([xStarNow(plotdims(1)),samplexNow(plotdims(1))]  ,   [xStarNow(plotdims(2)),samplexNow(plotdims(2))],'Color',colors(i,:));
%             end
%         end
%         if DVNow == maxDVValue
%             for k = 1:length(plotDimension)
%                 plotdims = plotDimension{k};
%                 max_dv_handle = plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'o','Color',[0,0,1]);
%             end
%         end
%         if DVNow == minDVValue
%             for k = 1:length(plotDimension)
%                 plotdims = plotDimension{k};
%                 min_dv_handle = plot(samplexNow(plotdims(1)),samplexNow(plotdims(2)),'o','Color',[1,0,1]);
%             end
%         end
%     end
%     drawnow
% end
% legend([max_dv_handle,min_dv_handle],sprintf('max DV = %f',maxDVValue), sprintf('min DV = %f',minDVValue))


















hold off

figure(rho_ROA_fig_Handle);
if SDSOS
    solverString = 'SDSOS Solver';
else
    solverString = 'SOS Solver';
end
suptitle(sprintf('ROA profile for %s Q=%s, stableEpsilon=%d, Qd factor=%d, PFL=%d',solverString,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,PFL))
% saveas(rho_ROA_fig_Handle,sprintf('%s/P_profile_Q_%s_ep_%d_Qd_factor_%d_SDSOS=%d.jpg',savePlotFolder,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,SDSOS)); 
saveas(rho_ROA_fig_Handle,sprintf('%s/P_profile_Q_%s_ep_%d_Qd_factor_%d_SDSOS=%d.pdf',savePlotFolder,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,SDSOS)); 
% saveas(rho_ROA_fig_Handle,sprintf('%s/P_profile_Q_%s_ep_%d_Qd_factor_%d_SDSOS=%d.fig',savePlotFolder,mat2str(diag(QFullNow)'),stableEpsilonNow,discreteQFactor,SDSOS)); 


end
end
end
end
end
