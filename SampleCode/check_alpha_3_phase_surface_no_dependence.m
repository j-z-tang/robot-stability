% This script plots the transversal surfaces obtained from the optimization
% in Step 3.

load alpha_3_abs_coord_unactuated/zResults_2015_08_04_8norm_normalMin_0.06_smoothCost_0_dz2dtCoef_0_dotProdCoef_0_pOpt_100.mat ...
    nominal_u_pp p Pipp PiDotpp zpp zDotpp xStarpp QFull f_full sampleStance stancePeriod stableEpsilon smoothCostCoef dz2dtCoef dotProdCoef

sampleNo = 40;
sampleTimeVec = 0:stancePeriod/(sampleNo-1):stancePeriod;

setAxis = true;
swingAxMin = -0.5;
swingAxMax = 0.5;
stanceAxMin = -0.45;
stanceAxMax = 0.45;

arrowSize = 0.1;
switchSurfSize = 0.2;
maxNormSurfSize = 1.5;
for i=1:sampleNo
    sampleTime = sampleTimeVec(i);
    xStarNow = xStarpp(sampleTime);
    uStarNow = nominal_u_pp(sampleTime);
    qSwAll(i,1)=xStarNow(1);
    qStAll(i,1)=xStarNow(2);
    qSwDotAll(i,1)=xStarNow(3);
    qStDotAll(i,1)=xStarNow(4);
    
    xStarAll(:,i) = xStarNow;
    uStarAll(:,i) = uStarNow;
    zppAll(:,i) = zpp(sampleTime); % each z a column
    zDotppAll(:,i) = zDotpp(sampleTime);
    fAll(:,i) = f_full(xStarNow,uStarNow);
end

% set up figure
figureHandle = figure();
f99Pos = get(figureHandle,'Position');
f99Pos = [f99Pos(1), f99Pos(2)-240, 880, 410];
set(figureHandle, 'Position', f99Pos)
orient(figureHandle,'landscape')
hold on;
zDirHandle=[];
colourMap = repmat([0,0.5,0.9],[sampleNo,1]);
for subplotIter = 1:2
figure(figureHandle); subplot(1,2,subplotIter); hold on
%% fill area of 'undefined dynamics' for original phase variable graph
startSwing = xStarAll(1,1); startStance = xStarAll(2,1);
endSwing = xStarAll(1,end); endStance=xStarAll(2,end);
if subplotIter==1
    % end point
    fill([endSwing,swingAxMax,swingAxMax,endSwing + (endStance-stanceAxMin)],[endStance,endStance,stanceAxMin, stanceAxMin],'y','EdgeColor','y')
    % beginning point
    fill([startSwing,swingAxMin,swingAxMin],[startStance,startStance, startStance + (startSwing-swingAxMin) ],'y','EdgeColor','y')
end


linesAll = [];
for i = 1:sampleNo
    t = sampleTimeVec(i);
    Pi = Pipp(t);
    x0 = xStarpp(t);
%     P = PMatrixpp(t);
    zNow = zppAll(:,i);
    zDotNow = zDotppAll(:,i);
    fNow = fAll(:,i);
    
    singPtDist = norm(zNow'*fNow)/norm(zDotNow);
    
    
    if singPtDist > maxNormSurfSize
        bigSingPtDist = true;
        normSurfSize = maxNormSurfSize;
    else
        bigSingPtDist = false;
        normSurfSize = singPtDist;
    end
    
    if subplotIter==1; zNow=[0,-1,0,0]'; normSurfSize = 1; end
    
    

    
    %% Plot direction of z
    zArrowBase = x0;
    zArrowTip = x0+(zNow*arrowSize);
    %     zDirHandle = line([zArrowBase(1),zArrowTip(1)],[zArrowBase(2),zArrowTip(2)],'Color',[1,0,0]);
    
    zDirHandle(i) = arrow(zArrowBase(1:2), zArrowTip(1:2));
    
    %% plot the normal
    % rotate tip 90 degrees clockwise
    normalDirClockwise = normSurfSize*[0, 1; -1, 0]*zNow(1:2);
    
    % rotate tip 90 degrees
    normalDirAntiClockwise = normSurfSize*[0, -1; 1, 0]*zNow(1:2);
    
    % sanity check to ensure normal
    if (~(zNow'*[normalDirClockwise;0;0] < 1e-10 )) || (~(zNow'*[normalDirAntiClockwise;0;0] < 1e-10))
        error('z not normal!')
    end
    
    clockwiseEdge = (x0(1:2)+normalDirClockwise)';
    antiClockwiseEdge = (x0(1:2)+normalDirAntiClockwise)';
    
    transSurfHandle = line([x0(1)+normalDirClockwise(1),x0(1)+normalDirAntiClockwise(1)],[x0(2)+normalDirClockwise(2),x0(2)+normalDirAntiClockwise(2)],'Color',colourMap(i,:));
    %% plot switching surface if applicable
    if i==1 || i==sampleNo
        switchSurfDirClockwise = switchSurfSize*[0, 1; -1, 0]*[-1;-1];
        switchSurfDirAntiClockwise = switchSurfSize*[0, -1; 1, 0]*[-1;-1];
        switchSurfHandle = line([x0(1)+switchSurfDirClockwise(1),x0(1)+switchSurfDirAntiClockwise(1)],[x0(2)+switchSurfDirClockwise(2),x0(2)+switchSurfDirAntiClockwise(2)],'Color','k','LineStyle','-.','LineWidth',1.1);
    end
    
end
invisibleRedLineforLegend=line([0.1,0.1],[0.1,0.1],'Color','r');

%% plot nominal trajectory
trajectoryHandle = plot(qSwAll,qStAll,'LineWidth',1.5,'color','k');
plot(qSwAll,qStAll,'x')

startptHandle = plot(qSwAll(1),qStAll(1),'og','LineWidth',3);
xlabel('Swing leg angle $\theta_{sw}$ (rads)','Interpreter', 'Latex')
ylabel('Stance leg angle $\theta_{st}$ (rads)','Interpreter', 'Latex')

axis equal;
view(90,-90) % swap the x and y axis
if setAxis
    axis([swingAxMin swingAxMax stanceAxMin stanceAxMax]);
end
if subplotIter==2
    title('New geometric phase variable')
else
    title('Traditional stance phase variable')
end
legend([switchSurfHandle,startptHandle,trajectoryHandle,invisibleRedLineforLegend,transSurfHandle],'Switching Surfaces','Starting Point','Nominal traj.','Phase Advancement Direction','Transversal Surface','Location','best')
% plot direction arrows for flow
arrow(zDirHandle,'EdgeColor','r','FaceColor','r','BaseAngle',60,'length',5);
noDirectionArrows = 2;
for arrowIter = 1:noDirectionArrows
    xStarIndex = floor(arrowIter*sampleNo/(noDirectionArrows+2));
    arrow(xStarAll([1,2],xStarIndex),xStarAll([1,2],xStarIndex+1));
end
end


clear figureHandle