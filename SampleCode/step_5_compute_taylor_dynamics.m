% This file computes the taylor series expansion of the dynamics at every
% sample point to enable SOS verification.  The result from this step is
% passed onto the SOS solver in the next step.

QFullList={diag([100 1 10 0.1])};

loadFolderRiccati = 'alpha_4_lyap_abs_coord_unactuated';
saveFolderTaylor = 'alpha_5_taylor_dynamics_abs_unactuated';

if ~exist(saveFolderTaylor,'dir')
    mkdir(saveFolderTaylor)
end

taylorOrder = 4;

for normalDeltaNow=[0.06]
    for smoothCostCoefNow=[0]
        for dz2dtCoefNow=[0]
            for dotProdCoefNow=[0]
                for stableEpsilonNow =  [1];
                    for discreteQFactor = [1]
                        for QFullIndex=1:length(QFullList)
                            QFullNow=QFullList{QFullIndex};
                            pOpt = 100;
                            clock = tic;
                            
                            % check if done before...
                            if exist(sprintf('%s/spline_full_dynamics_taylor_stableEp_%d_Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d.mat',saveFolderTaylor,stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor),'file')
                                prompt = 'Previous computation found.  Recompute Taylor expansions?  Note if the .mat files were moved from another machine, you must recompute them.  Y/N [N]';
                                userInput = input(prompt,'s');
                                if strcmp(userInput,'Y')
                                    display('Recomputing Taylor expansions...')
                                else
                                    return
                                end
                            end
                            
                            % try loading data from previous steps
                            try
                                % load z_Pi_P_splines_extra_stable_ep_1__Q_100_50_10_5__ZnormalMin_0.10_ZsmoothCost_10_dz2dtCoef_8.mat PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull stateAll sampleStance limCycDuration dfduNowpp normalDelta smoothCostCoef dz2dtCoef
                                eval(strrep(sprintf('load  %s/z_Pi_P_bezierTheta_extra_stable_ep_%d__Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt%d_Qd_factor_%d.mat p nominal_u PMatrixpp PDotMatrixpp Pipp PiDotpp zpp zDotpp xStarpp f QFull sampleStance limCycDuration dfduNowpp normalDelta smoothCostCoef dz2dtCoef symImpactMap dotProdCoef stableEpsilon',...
                                    loadFolderRiccati,stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,pOpt,discreteQFactor),'/',filesep))
                            catch err
                                display(sprintf('COULD NOT LOAD  pre_2_1_mat/z_Pi_P_splines_extra_stable_ep_%d__Q_%s_biggerImpactQ__ZnormalMin_%.2f_ZsmoothCost_%d_dz2dtCoef_%d_dotProdCoef_%d.mat',...
                                    stableEpsilonNow,mat2str(diag(QFullNow)'),normalDeltaNow,smoothCostCoefNow,dz2dtCoefNow,dotProdCoefNow))
                                continue
                            end
                            
                            
                            
                            QTrans = @(t) Pipp(t)*QFull*(Pipp(t)');
                            RInverse = 100;
                            sampleNo = 40;
                            sampleTimeVec = 0:limCycDuration/(sampleNo-1):limCycDuration;
                            fSymbolicEvalWithInput = f;
                            stateAll = sampleStance.q;
                            
                            %% construct discrete transverse dynamics
                            syms qSw qSt qSwDot qStDot
                            impactPt = stateAll(:,end);
                            startPt = stateAll(:,1);
                            xper = sym('xper', [3 1]); transState = [xper(1); xper(2); xper(3)];
                            
                            
                            %% construct continuous FULL dynamics
                            uCon = sym('uCon','real');
                            x=sym('x',[4 1]);
                            x = sym(x,'real');
                            fFullTaylorSample = cell([sampleNo, 1]);
                            fFullTaylorInxperpSample = cell([sampleNo, 1]);
                            xperDotOriginal = cell([sampleNo, 1]);
                            xperDotTaylorAll = cell([sampleNo, 1]);
                            xperDotOriginalDiff = cell([sampleNo, 1]);
                            dfperdxperTaylor = cell([sampleNo, 1]);
                            QTransAll= cell([sampleNo, 1]);
                            BtransAll= cell([sampleNo, 1]);
                            PMatrixAll= cell([sampleNo, 1]);
                            PBRInverseBP= cell([sampleNo, 1]);
                            
                            %% switching condition taylor
                            xStarEnd = xStarpp(sampleTimeVec(end));
                            [deltaq, deltaqDot] = symImpactMap(x,p);
                            deltaqDotTaylor = taylor(deltaqDot, x, xStarEnd,'Order',taylorOrder);
                            
                            %% transverse dynamics taylor
                            % for i = 1:sampleNo
                            parfor i = 1:sampleNo
                                %if ~mod(i,5)
                                display(sprintf('Calculating taylor series expansion of dynamics for transversal surface no. %d', i))
                                %end
                                
                                timeNow = sampleTimeVec(i);
                                
                                xStarNow = xStarpp(timeNow);
                                uStarNow = nominal_u(timeNow);
                                PiNow = Pipp(timeNow);
                                PiDotNow = PiDotpp(timeNow);
                                zNow = zpp(timeNow);
                                zdotNow = zDotpp(timeNow);
                                
                                fNow = fSymbolicEvalWithInput(x,uCon);
                                fNowTaylor = taylor(fNow,[x;uCon],[xStarNow;uStarNow],'Order',taylorOrder);  % With uStar here.
                                fFullTaylorSample{i}=fNowTaylor;
                                
                                xInxperp = xStarNow+PiNow.'*xper;
                                fNowInxperp = fSymbolicEvalWithInput(xStarNow+PiNow.'*xper, uCon);
                                
                                % transverse control
                                PMatrixNow = PMatrixpp(timeNow);
                                dfduNow = dfduNowpp(xStarNow(1),xStarNow(2),xStarNow(3),xStarNow(4));
                                fEvalxStar = fSymbolicEvalWithInput(xStarNow,uStarNow); % with ustar here
                                dtaudu = zNow'*dfduNow/(zNow'*fEvalxStar);
                                transverseB = PiNow*dfduNow - PiNow * fEvalxStar * dtaudu;
                                
                                u = uStarNow + (-1)*RInverse*transverseB'*PMatrixNow*xper;
                                
                                fNowInxperpWithFeedback = subs(fNowInxperp, uCon, u);
                                fNowInxperpWithFeedbackTaylor = taylor(fNowInxperpWithFeedback, transState, 'Order', taylorOrder);
                                fFullTaylorInxperpSample{i} = fNowInxperpWithFeedbackTaylor;
                                
                                taudot = (zNow'*fNowInxperpWithFeedback)/(zNow'*fEvalxStar - zdotNow'*PiNow'*xper);
                                xperdot = taudot*PiDotNow*PiNow'*xper + PiNow*fNowInxperpWithFeedback ...
                                    - PiNow*fEvalxStar*taudot;
                                
                                xperDotOriginal{i} = xperdot;
                                xperDotOriginalDiff{i} = jacobian(xperdot);
                                dfperdxperTaylor{i} = taylor( xperDotOriginalDiff{i} , transState, 'Order', taylorOrder);
                                xperDotTerm1Taylor=taylor(taudot*PiDotNow*PiNow'*xper,transState,'Order',taylorOrder);
                                xperDotTerm2Taylor=taylor(PiNow*fNowInxperpWithFeedback,transState,'Order',taylorOrder);
                                xperDotTerm3Taylor=taylor(-(PiNow*fEvalxStar*taudot),transState,'Order',taylorOrder);
                                
                                xperDotTaylorAll{i} = xperDotTerm1Taylor+xperDotTerm2Taylor+xperDotTerm3Taylor;
                                
                                QTransAll{i} = QTrans(timeNow);
                                BtransAll{i} = transverseB;
                                PMatrixAll{i} = PMatrixNow;
                                PBRInverseBP{i} = PMatrixNow * transverseB * RInverse * transverseB' * PMatrixNow;
                            end
                            
                            
                            
                            finished_taylor_full_dynamics = toc(clock)
                            eval(sprintf('save  %s/spline_full_dynamics_taylor_stableEp_%d_Q_%s_biggerImpactQ__ZnormalMin_%.2f_pOpt_%d_Qd_factor_%d.mat',saveFolderTaylor,stableEpsilonNow,mat2str(diag(QFull)'),normalDeltaNow,pOpt,discreteQFactor))
                        end
                    end
                end
            end
        end
    end
end
