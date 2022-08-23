function [ x2 ] = DynamicAlignDev( y1, y2, frameWidth, frameStep)
%DYNAMICALIGNDEV Aligns a signal y2 to signal y1 by cross-correlating frames of y2 which are frameWidth samples wide, across all of y2
% frameStep controls how far the frame is moved each time (e.g. frameStep==frameWidth means there is no overlap). I tend to use frameStep = 0.5*frameWidth

%THIS IS THE DEV VERSION - NOT OPTIMISED FOR PERFORMANCE

tic;

len1 = ceil((length(y1) - frameWidth));
len2 = ceil(length(y2) / frameStep);

npos = ceil(length(y1) ./ frameStep) - 1;

figure(201); hold off;
plot(y1);
title('y1');

% figure(200); hold off;

CM = zeros(length(y1), npos);
P = zeros(npos, 1);

for i=1:npos
    
    %p1 = ((i-1)*frameStep) + 1:frameWidth; 
    p2 = ((i-1)*frameStep) + (1:frameWidth);
    [~, ii] = find(p2 > length(y2), 1);
    if ~isempty(ii)
        if ii == 1
            p2 = [];
        else
            p2 =  p2(1:ii-1);
        end
    end
        
    if ~isempty(p2)
        [C, lags] = xcorr(y1, y2(p2));  %Compare window of y2 across all of y1
        %Result length [length(C)] is == length(y1)*2-1

        ci = (length(y1)+1+(i-1)*frameStep):length(C);
        if length(ci) == length(y1)
            CM(:,i) = C(ci);
        else
            ovr = length(y1) - length(ci);
            CM(:,i) = [C(ci) ; zeros(ovr, 1)]; %Pad overrun with zeros
        end 

        [~, mp] = max(C(ci));
        P(i) = mp;

        %CM(:,i) = C;

    %     figure(202);
    %     hold off;
    %     plot(y2(p2));
    %     title(['y2 window, i=' num2str(i)]);
    % 
    %     figure(200);
    %     plot(lags, C);
    %     hold on;
    %     title(['C score, i=' num2str(i)]);

    end
end
toc;

%% --- Find the best path from one side to the other ---

tic;
score = zeros(length(y1), 1);
routes = zeros(length(y1), npos);
for start=1:length(y1)
    
    score(start) = CM(start, 1);
    lag = start;
    routes(start, 1) = start;
    
    for p=2:npos    
        
        s0 = CM(lag, p); %Favour central path
        if lag > 1
            sL = CM(lag-1, p) ./ 1.414;
        else
            sL = -inf;
        end
        if lag < length(y1)
            sR = CM(lag+1, p) ./ 1.414;
        else
            sR = -inf;
        end
        
        if (sL > sR) && (sL > s0)
            sB = sL;
            lag = lag - 1;
        elseif sR > s0
            sB = sR;
            lag = lag + 1;
        else
            sB = s0;
        end;
        
        score(start) = score(start) + sB;
        routes(start, p) = lag;
    end
end

[~, bestScoreIndex] = max(score)

toc;


%% --- Plot dat shi' ---

figure(204); hold off;
surf(0:length(y1)-1, 0:npos-1, CM', 'LineStyle', 'none');
view(0, 90);
xlabel('lag');
ylabel('window position');

colormap(genCM('BlueWhiteRed'));
c=colorbar;

hold on;
plot3(P, 0.5 + (0:npos-1), ones(npos, 1) * 1, 'k+');
hold on;
plot3(routes(bestScoreIndex,:), 0.5 + (0:npos-1), ones(npos, 1) * 1, 'k');