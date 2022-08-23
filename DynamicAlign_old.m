function [ x2shift, x2score, lagScores ] = DynamicAlign( y1, y2, frameWidth, frameStep, alignRange, plotGraphs)
%DYNAMICALIGNXC Aligns a signal y2 to signal y1 by calculating the squares of differences between frames of y2, across all of y1
% frameWidth is the width of each frame of y2, in samples
% frameStep controls how far the frame is moved each time (e.g. frameStep==frameWidth means there is no overlap). I tend to use frameStep = 0.5*frameWidth
% alignRange sets how many samples either side of the zero position the frame should be compared at - i.e. the maximum x difference between signals y1 and y2
% plotGraphs (optional, default 0) set to 1 to plot a graph of the alignment output
% RETURNS x2shift contains values corresponding to each y2 value, indicating how far to shift each sample on the x axis (in non-integer samples)
% RETURNS x2score contains confidence scores pertaining to each x2shift value
% RETURNS lagScores contains overall scores at each lag position
%
% If you want to return a sample-aligned version of y2, use DynamicAlignResamp instead

if ~exist('plotGraphs')
    plotGraphs = 0;
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
clen = alignRange*2+1; %Comparison length

CM = zeros(clen, npos); %score matrix
P = zeros(npos, 1);

%### THIS MAY BENEFIT FROM PARFOR IN THE FUTURE BUT I HAVEN'T GOT AROUND TO OPTIMISING IT YET
%tic;
%ticBytes(gcp);
for i=1:npos
    
    fpos = (i-1)*frameStep; %Frame position of element 1 : Single position in y2, Central position (of k) in y1
    
    % Sample indices for y1 and y2 data
    %   in the case of y1 this is at the central position (i.e. when k == 0)
    p = fpos + (1:frameWidth);
    vpi = find((p <= L1) & (p <= L2), 1, 'last'); %Find last valid index of p (not out of range for L1 or L2)
    p = p(1:vpi);
    
    %Cross correlate manually so we can control range
    C = zeros(clen, 1);
    for k = -alignRange:alignRange
        c = 0;
        p1_k = p + k; %Get shifted indices for y1 dataset

        for m = 1:length(p)
            if (p1_k(m) > 0) & (p1_k(m) <= L1)
                %c = c + y1(p1_k(m)) .* y2(p(m)); %Cross-correlation
                %c = c + 1./(abs(y1(p1_k(m)) - y2(p(m)))+1).^2; %Squares of differences
                c = c + (y1(p1_k(m)) - y2(p(m))).^2; %Squares of differences
            end
        end
        
        C(k+alignRange+1) = c;
    end    
    CM(:,i) = C;
    
    [~, mp] = max(C);
    P(i) = mp;
end
%tocBytes(gcp);
%toc;

%% --- Find the best path from one side to the other ---

%tic;
score = zeros(clen, 1);
routes = zeros(clen, npos);
routescores = zeros(clen, npos);
for start=1:clen
    
    score(start) = CM(start, 1);
    lag = start;
    routes(start, 1) = start - alignRange - 1;
    
    %For each step we can step left, right, or straight forward
    % We calculate scores for these, sL, sR, and s0, and choose the best
    for p=2:npos    
        
        s0 = CM(lag, p); %Favour central path
        if lag > 1
            sL = CM(lag-1, p) ./ 1.414;
        else
            sL = NaN;
        end
        if lag < clen
            sR = CM(lag+1, p) ./ 1.414;
        else
            sR = NaN;
        end
        
        if (sL < sR) && (sL < s0) && (~isnan(sL)) && lag > 1
            sB = sL;
            lag = lag - 1;
        elseif sR < s0 && (~isnan(sR)) && lag < clen
            sB = sR;
            lag = lag + 1;
        else
            sB = s0;
        end;
        
        score(start) = score(start) + sB;
        routes(start, p) = lag - alignRange - 1; %This gets the lag in range -alignRange to +alignRange (rather than 1:clen as lag is currently)
        routescores(start, p) = sB; %Keep track of score for each step of each route - for use in calculating the error
    end
end

[~, bestScoreIndex] = min(score);

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
    surf(-alignRange:alignRange, 0:npsurf-1, CM', 'LineStyle', 'none');
    view(0, 90);
    xlabel('lag');
    ylabel('window position');

    addpath('../IMUMatlab');
    colormap(genCM('BlueWhiteRed'));
    colorbar;

    hold on;
    plot3(P-alignRange, 0.5 + (0:npos-1), ones(npos, 1) * max(CM) + 1, 'k+');
    hold on;
    plot3(routes(bestScoreIndex,:), 0.5 + (0:npos-1), ones(npos, 1) * max(CM) + 1, 'k');
    
    disp(['Best Score Index = ' num2str(bestScoreIndex)]);
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

