function [ reachStart, reachStop, reachMW, pos1, pos2 ] = getReachTimesJPBeth1st( positionData )
%Extracts best reach start times from .ns2 data unless BabData set to 1 for
%processed ContData files.
% Can be made to extract from other data types. Assumes 1kHz sampling
%   Assumes  X and Y position are channels 2 and 3 (line 13)

%outputs:
%  reachStart and reachStop are vectors of reach times
%  reachMW is the amplitude readout of the whole session
%  pos1 = all reach traces, aligned to start
%  pos2 = all reach traces, aligned to stop
% 
%  sgolay injects a negative dip before a reach as an artifact of the filter.
%  we can use that for reach detection, but makes for poor position traces
%  sgolay filter is a convolution filter which fits successive sub-sets of
%  adjacent data points with a low-degree polynomial by the method of
%  linear least squares. 

%TbetweenR = 2000; % minimum time between reaches, in ms
TbetweenR = 3000; % minimum time between reaches, in ms

maxDur=1500;      % max reach duration

NsDATA.Data(2,:) = positionData(1,:); % Xposition (Beth = trial starts)
NsDATA.Data(3,:) = positionData(2,:); % Yposition (Beth = wheel)

reachSG = sgolayfilt(sqrt((NsDATA.Data(3,:)-median(NsDATA.Data(3,:))).^2 ... 
    + (NsDATA.Data(2,:)-median(NsDATA.Data(2,:))).^2),3,201); % smoothing with sgolay filter
reachMW = moveavg(sqrt((NsDATA.Data(3,:)-median(NsDATA.Data(3,:))).^2 ...
    + (NsDATA.Data(2,:)-median(NsDATA.Data(2,:))).^2),10);    % smoothing with moving average filter

vel=[0, diff(reachSG)]; % velocity - differentiation of the position trajectory

scc=sort(reachSG);   % sort position data
scd=sort(vel);      % sort velocity data 
scc(scc<median(scc))=median(scc);
cutoff  = scc(round(length(scc)*.90));  % set the 90% cutoff for the position data
cutoffd = scd(round(length(scd)*.98));  % set the 98% cutoff for the velocity data (velocity cutoff is required, since reaches are expected to occur at high velocity)
allValidSamps1=find(vel>cutoffd);       % velocities exceeding the cutoff
allValidSamps1=allValidSamps1(allValidSamps1<(length(reachSG)-300)); % this line is just to exclude the data points at the tail
allValidSamps=[];

for jj = 1:length(allValidSamps1) % all points exceeding the velocity cutoff
    if length(allValidSamps)==0   % in case of the 1st high-velocity data point
        if sum(reachSG(allValidSamps1(jj):allValidSamps1(jj)+200)>cutoff)>20 % this part is just to sort out high-velocity reaches from the short jerks which is not really a reach
           allValidSamps=[allValidSamps, allValidSamps1(jj)];               % if it satisfies the condition for a legitimate reach, then register it as a reach
        end;
    else    % from the 2nd high-velocity data point on 
        if (allValidSamps1(jj)-allValidSamps1(jj-1))>100           % if the current high-velocity point came 100 ms or more after the previous high-velocity point; this will prevent the redundant counting from a continued reach
            if (allValidSamps1(jj)-allValidSamps(end)) > TbetweenR % this is to ensure that the interval between reaches is greater than the set interval threshold (e.g. 2000 ms)
                if sum(reachSG(allValidSamps1(jj):allValidSamps1(jj)+200)>cutoff)>20 % finally ensure that it is a well-timed reach not a short jerk
                    allValidSamps=[allValidSamps, allValidSamps1(jj)];              % if it satisfies conditions for a legitimate reach, then register it as a reach
                end
            end
        end
    end
end

%% detect reach start
reachStart=[]; % the goal here is to find the near-zero point that is immediately prior to each high-velocity point  
for i  = 1:length(allValidSamps)
    ii = allValidSamps(i);                      % current valid sample (high-velocity data point) 
    iii= find(reachSG(1:ii)<cutoff*.4,1,'last'); % find the most recent near zero-point
    if iii>10
        if ii-iii>1000  % sometimes the joystick doesn't reset to zero and induces 20 second long reaches. this prevents that by setting the second (additional) cutoff
            scc2=sort((reachSG(iii:ii)));
            cutoff2 = scc2(round(length(scc2)*.90)); % reset the cutoff2 to be a point immediatly before the reach peak (farthest point)
            iii=find(reachSG(1:ii)<cutoff2,1,'last');
        end
        realStartMinus10=find(vel(1:iii)<0,10,'last'); % find 10 elements whose velocity is negative prior to the current near zero point  
    else % if there's no near zero point
        realStartMinus10=ii; % take the current valid sample (high-velocity data point) as a reach start point
    end
    reachStart = [reachStart, realStartMinus10(1)];
end

reachStart=reachStart([1,find(diff(reachStart)>TbetweenR)+1]); % ensure that intervals between reaches are greater than 2 sec
reachStart=unique(reachStart)+10; % unique ensures that there's no redundant reaches

% this will eliminate pseudo starts, you might not want to do eliminate them
pos1=[];
for i = reachStart
    if i > 1000
        if i+1499 > length(reachMW) % in case the reach trajectory end exceeds the length of the time series (this must happen only for the last reach)
            %pos1=[pos1; reachMW(i-200:end)-reachMW(i)]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
        else % in case the reach trajectory is within the length of the time series
            pos1=[pos1; reachMW(i-200:i+1499)-reachMW(i)]; % position data corresponding to each reach from 200-ms before to 1500-ms after each reach start
        end
    end
end
% figure; imagesc(pos1); colorbar
% title('start times, aligned at 200ms')
      
%% detect reach stop
reachStop=[]; % the goal is to find the reach below the threshold for a certain duration of time (150 ms)
for i = 1:length(reachStart)
    
    if reachStart(i)+maxDur < length(reachMW)
    ii= reachMW(reachStart(i):reachStart(i)+maxDur); % take the smoothed position data corresponding to each reach 
    end
    
    [~,reachPeakTime]=max(ii(1:700));   % detect reach peak (the farthest reach point)
    ii2=ii(reachPeakTime:end);          % peak-to-end of reach  
    sortStop=sort(ii2);
    binSize=5;  % 10 may be better
    jsBinned=hist(ii2,binSize); %jsBinned=hist(ii2,binSize); % jsBinned=hist(ii2,min(ii2):binSize:max(ii2));
    [~,cutoffStop1]=max(jsBinned(1:round(length(jsBinned)/3)));
    cutoffStop = sortStop(jsBinned(cutoffStop1)); % min(ii2) + sortStop(jsBinned(cutoffStop1)); % the way I think the cutoffStop should be defined!
    %cutoffStop=cutoffStop1*binSize+10+min(ii2); % cutoffStop1
    
    N=150; % # of ms consecutively below threshold needed to detect a stable reach stop point
    t = find(ii2<cutoffStop);  % find the data points below the cutoffStop 
    x = diff(t)==1;            % spot out the points consecutively below the cutoffStop
    f = find([false,x]~=[x,false]); % spot out the non-consecutive data points 
    g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first'); % find the first data point, at which the interval between non-consecutive points is greater than 150 ms, as this means that there are 150 or more points consecutively under the reachstop threshold
    almostEnd=t(f(2*g-1))-1;   % just tranlate g to a point on t  
    if length(almostEnd)==0    % in case there's no such point
        almostEnd=find(ii2<cutoffStop,1);   % just take the point below the reach threshold
        %reachStart(i);
    end
    %thisStop= find(moveavg(diff(ii2(almostEnd:end)),10)>=0,1);
    reachStop=[reachStop, reachPeakTime+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
end

pos2=[];
for i =reachStop
    if i>1500 && i < (length(reachMW)-200)
        pos2=[pos2; reachMW(i-1499:i+200)-reachMW(i)];
    end
end

% figure; imagesc(pos2); colorbar
% title('stop times, aligned at 1500ms')

% figure; plot(reachMW);
% hold on; plot(reachStart, reachMW(reachStart),'g*');
% hold on; plot(reachStop, reachMW(reachStop),'r*');
% title('JS position with start and stop times');

% Trial-by-trial inspection 
% reachNumb = 152; % trial
% figure; plot(pos1(reachNumb,:)); hold on; 
% plot(200, 0, 'g*'); % the reach start point is 200th on the pos1, and set to 0 (line 82)
% plot(200 + reachStop(reachNumb) - reachStart(reachNumb), reachMW(reachStop(reachNumb)) - reachMW(reachStart(reachNumb)), 'r*'); % the reach start point is 200th on the pos1, and set to 0 (line 82)

%%
      %  if max(diff(endThresh))==1
           % almostEnd=endThresh(end)+reachStart(i);
            
       % else
        %iicut=find(diff(find(reachMW>cutoffStop))>
        %end
       % reachStop=[reachStop,find(diff(reachMW(almostEnd:...
        %    almostEnd+maxDur))>=0,1)+almostEnd];
end
   % below displacement threshold and 