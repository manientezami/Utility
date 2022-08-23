function [ x2shift, x2score, lagScores ] = DynamicAlign( y1, y2, frameWidth, frameStep, alignRange, plotGraphs, crossFrameLagChange, zeroLagIndex)
%DYNAMICALIGNXC Aligns a signal y2 to signal y1 by calculating the squares of differences between frames of y2, across all of y1
% frameWidth is the width of each frame of y2, in samples
% frameStep controls how far the frame is moved each time (e.g. frameStep==frameWidth means there is no overlap). I tend to use frameStep = 0.5*frameWidth
% alignRange if 2-elements sets minimum and maximum samples either side of the zero position the frame should be compared at - i.e. the maximum x difference between signals y1 and y2 - if 1-element the range is -alignRange:alignRange
% plotGraphs (optional, default 0) set to 1 to plot a graph of the alignment output
% crossFrameLagChange (optional, default 1) sets the max allowed lag change from frame position to frame position
% zeroLagIndex (optional, default 1) controls where in y1 the alignRange zero position is
% RETURNS x2shift contains values corresponding to each y2 value, indicating how far to shift each sample on the x axis (in non-integer samples)
% RETURNS x2score contains confidence scores pertaining to each x2shift value
% RETURNS lagScores contains overall scores at each lag position
%
% If you want to return a sample-aligned version of y2, use DynamicAlignResamp instead

if ~exist('plotGraphs')
    plotGraphs = 0;
end
if ~exist('crossFrameLagChange')
    crossFrameLagChange = 1;
end
if ~exist('zeroLagIndex')
    zeroLagIndex = 1;
end

if length(y2) < frameStep
    frameWidth = length(y2);
end

L1 = length(y1);
L2 = length(y2);
npos = ceil(length(y1) ./ frameStep) - 1; %Number of possible frame positions
if npos < 1
    npos = 1; %### This feels bodgey -- not sure
end

if length(alignRange) == 2
    rngN = alignRange(1);
    rngP = alignRange(2);
    if rngP < rngN
        error('alignRange second element should be higher than the first');
    end
elseif length(alignRange) == 1
    rngN = -alignRange;
    rngP = alignRange;
else
    error('alignRange should have 1 or 2 elements');
end
    
clen = rngP-rngN+1; %Comparison length

CM = zeros(clen, npos); %score matrix
P = zeros(npos, 1);

%### THIS MAY BENEFIT FROM PARFOR IN THE FUTURE BUT I HAVEN'T GOT AROUND TO OPTIMISING IT YET
%tic;
%ticBytes(gcp);
C = zeros(clen, 1);
for i=1:npos
    
    fpos = (i-1)*frameStep; %Frame position of element 1 : Single position in y2, Central position (of k) in y1
    
    % Sample indices for y1 and y2 data
    %   in the case of y1 this is at the central position (i.e. when k == 0)
%     p = fpos + (1:frameWidth);
%     vpi = find((p <= L1) & (p <= L2), 1, 'last'); %Find last valid index of p (not out of range for L1 or L2)
%     p = p(1:vpi);
%     LP = length(p);
    
    %Cross correlate manually so we can control range
    %C = zeros(clen, 1);
    %cv = zeros(length(p), 1);
    y2start0 = fpos+1;
    y2end0 = min(y2start0+frameWidth-1, L2);
    for k = rngN:rngP
        %c = 0;
        
        y2start = y2start0;
        y2end = y2end0;
        y1start = y2start + k + (zeroLagIndex-1);
        y1end = y2end + k + (zeroLagIndex-1);
        %y1start = p1_1;
        %y1end = min(y1start+L1-1, L1);
        %y2end = min(y2start+L2-1, L2);
        %p2start = 2-p1_1;
        %p2end = min(p2start+L1-1, LP);
        
%         if p2start < 1
%             %y1start = y1start - p2start + 2;
%             p2start = 1;
%         end
%         if y1start < 1
%             %p2start = p2start - y1start + 2;
%             y1start = 1;
%         end

%         y1len = y1end-y1start+1;
%         p2len = p2end-p2start+1;
%         dlen = p2len - y1len;
        
        if y1start < 1
            y2start = y2start - y1start + 1;
            y1start = 1;
        end
        if y2start < 1
            y1start = y1start - y2start + 1;
            y2start = 1;
        end
        if y1end > L1
            y2end = y2end - (y1end-L1);
            y1end = L1;
        end
        if y2end > L2
            y1end = y1end - (y2end-L2);
            y2end = L2;
        end
        
%         y1len = y1end-y1start+1;
%         y2len = y2end-y2start+1;
%         dlen = y2len - y1len;
%         if dlen < 0 %y1len is longer
%             y1end = y1end + dlen;
%         elseif dlen > 0 %y2len is longer
%             y2end = y2end - dlen;
%         end        
% %         if dlen < 0 %y1len is longer
% %             y1end = y1end + dlen;
% %         elseif dlen > 0 %p2len is longer
% %             p2end = p2end - dlen;
% %         end
        
        %dylen = p2end-p2start+1;
        %p1start = p1_1;
        %p1end = p1_k(mend);
        
        %dy = y1(y1start:y1end) - y2(p(p2start:p2end));
        dy = y1(y1start:y1end) - y2(y2start:y2end);
        in = ~isnan(dy);
        c = sum(dy(in).^2);
        if length(in) > 0
            c = c ./ length(in); %"Normalise" scores
        end
        
%         p1_k = p + k + (zeroLagIndex-1); %Get shifted indices for y1 dataset
%         ctr = 0;
%         for m = 1:length(p)
%             if (p1_k(m) > 0) & (p1_k(m) <= L1)
%                 %c = c + y1(p1_k(m)) .* y2(p(m)); %Cross-correlation
%                 %c = c + 1./(abs(y1(p1_k(m)) - y2(p(m)))+1).^2; %Squares of differences
%                 if ~isnan(y1(p1_k(m)))
%                     if ~isnan(y2(p(m)))
%                         c = c + (y1(p1_k(m)) - y2(p(m))).^2; %Squares of differences    
%                     else %use y1 only to penalise
%                         c = c + y1(p1_k(m)).^2; %Squares of y1 (difference between y1 and 0)
%                     end
%                 else
%                     if ~isnan(y2(p(m))) %use y2 only to penalise
%                         c = c + y2(p(m)).^2; %Squares of y2 (difference between y2 and 0)
%                     %else
%                         %c = c + 0; %add difference between 0 and 0
%                     end
%                 end
%                 ctr = ctr+1;
%             else %use y2 only to penalise - We dont do this anymore - because it means with equal length data (y1 and y2 equal length), positive shifts are penalised more - the y1 overhangs (at the start) are not penalised
%                 %c = c + y2(p(m)).^2; %Squares of y2 (difference between y2 and 0)
%             end
%             %cv(m) = c;
%         end
        
%         if c ~= c2
%             disp(['!!! c not equal to c2, k=' num2str(k) ', i=' num2str(i)]);
%         end
% 
%         if mod(k, 1000) == 1
%             disp(['k = ' num2str(k) ' : ' num2str(c) ' , ' num2str(c2)]);
%         end
        
%         if k==5500 || k == 500
%             figure(k);
%             hold off;
%             plot(diff(cv));
%             title(['k=' num2str(k) ', score=' num2str(cv(end))]);
%         end
        
        %C(k+alignRange+1) = c;
        CM(k-rngN+1,i) = c;
    end    
    %figure(8); plot(-alignRange:alignRange, C)
    %CM(:,i) = C;
    
    %[~, mp] = max(C);
    [~, mp] = max(CM(:,i));
    P(i) = mp;
end
%tocBytes(gcp);
%toc;

%% --- Find the best path from one side to the other ---

%tic;
score = zeros(clen, 1);
routes = zeros(clen, npos);
routescores = zeros(clen, npos);
changeRange = -crossFrameLagChange:crossFrameLagChange;
for start=1:clen
    score(start) = CM(start, 1);
    lag = start;
    routes(start, 1) = start + rngN - 1;
    
    %For each step we can step left, right, or straight forward
    % We calculate scores for these, sL, sR, and s0, and choose the best
    for p=2:npos    
        
        localScores = NaN(1, crossFrameLagChange*2 + 1);
        
        for lc = -crossFrameLagChange:crossFrameLagChange
            localLag = lag + lc;
            if localLag >= 1 && localLag <= clen
                localScores(lc+crossFrameLagChange+1) = CM(localLag, p);
            end
        end
        
        lsMinVal = min(localScores);
        lsMinIdx = find(localScores == lsMinVal); %Doing it like this because there could be >1 minimum value
        if length(lsMinIdx) > 1 %More than one min value
            [~, lsmiBestIdx] = min(abs(changeRange(lsMinIdx))); %Find the min position with minimum change
            lsMinIdx = lsMinIdx(lsmiBestIdx); %Narrow down to one option
        end
        lag = lag + changeRange(lsMinIdx); %Set the new lag
            
        score(start) = score(start) + lsMinVal;
        routes(start, p) = lag + rngN - 1; %This gets the lag in range -alignRange to +alignRange (rather than 1:clen as lag is currently)
        routescores(start, p) = lsMinVal; %Keep track of score for each step of each route - for use in calculating the error
        
%         s0 = CM(lag, p); %Favour central path
%         if lag > 1
%             sL = CM(lag-1, p) ./ 1.414;
%         else
%             sL = NaN;
%         end
%         if lag < clen
%             sR = CM(lag+1, p) ./ 1.414;
%         else
%             sR = NaN;
%         end
%         
%         if (sL < sR) && (sL < s0) && (~isnan(sL)) && lag > 1
%             sB = sL;
%             lag = lag - 1;
%         elseif sR < s0 && (~isnan(sR)) && lag < clen
%             sB = sR;
%             lag = lag + 1;
%         else
%             sB = s0;
%         end;
%         
%         score(start) = score(start) + sB;
%         routes(start, p) = lag - alignRange - 1; %This gets the lag in range -alignRange to +alignRange (rather than 1:clen as lag is currently)
%         routescores(start, p) = sB; %Keep track of score for each step of each route - for use in calculating the error
    end
end

bestScore = min(score);
bestScoreIndex = find(score == bestScore); %Doing it like this because there could be >1 minimum value
if length(bestScoreIndex) > 1 %More than one min value
    arng = rngN:rngP;
    [~, bbsIdx] = min(abs(arng(bestScoreIndex))); %Find the min position with minimum realignment (lowest alignRange value)
    bestScoreIndex = bestScoreIndex(bbsIdx); %Narrow down to one option
end

%toc;



%% Plottage

if plotGraphs == 1
    npsurf = npos;
    if size(CM,2) == 1
        CM = [CM, CM]; %Bodge for display only
        npsurf = npos + 1;
    end
    
    %figure(204); 
    hold off;
    surf(rngN:rngP, 0:npsurf-1, CM', 'LineStyle', 'none');
    view(0, 90);
    xlabel('lag');
    ylabel('window position');

    addpath('../IMUMatlab');
    colormap(genCM('BlueWhiteRed'));
    colorbar;

    hold on;
    plot3(P+rngN, 0.5 + (0:npos-1), ones(npos, 1) * max(CM) + 1, 'k+');
    hold on;
    plot3(routes(bestScoreIndex,:), 0.5 + (0:npos-1), ones(npos, 1) * max(CM) + 1, 'k');
    
    %disp(['Best Score Index = ' num2str(bestScoreIndex)]);
    arng = rngN:rngP;
    disp(['Best score lag = ' num2str(arng(bestScoreIndex))]);
end


%% --- Calculate X values ---

xin = (0:npos-1)*frameStep;
xout = 0:L2-1;

if length(xin) > 1
    x2shift = interp1(xin, routes(bestScoreIndex,:), xout, 'pchip')';
    x2score = interp1(xin, routescores(bestScoreIndex,:), xout, 'pchip')';
else
    x2shift = ones(length(xout), 1) .* routes(bestScoreIndex,:);
    x2score = ones(length(xout), 1) .* routescores(bestScoreIndex,:);
end

%x2shift = movmean(x2shift, frameWidth*2+1); %Smooth it out a bit
% wl = frameWidth*16;
% [b,a] = butter(3, 2/wl, 'low');
% x2shift = filtfilt(b,a, x2shift); %Smooth it out a bit

%x2shift = interp1(xin, routes(bestScoreIndex,:), xout);
%x2score = interp1(xin, routescores(bestScoreIndex,:), xout);
if xout(end) > xin(end)
    %x2shift(xin(end):end) = x2shift(xin(end)); %Hold end value
    %x2score(xin(end):end) = x2score(xin(end)); %Hold end value
    
    x2shift(xin(end)+1:end) = x2shift(xin(end)+1); %Hold end value
    x2score(xin(end)+1:end) = 0; %Set to zero
end

if plotGraphs == 1
    plot3(x2shift, (0:L2-1) / frameStep + 0.5, ones(L2, 1) .* max(CM) + 1, 'm');
end

lagScores = score;

