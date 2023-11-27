%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To check the duplicates (MUs that share most of their discharge times,
% the threshold is defined in the script main.m)

% Input: 
%   PulseT = Pulse train of each MU
%   distime = discharge times of the motor units 
%   muscle = indexes of the motor units 
%   maxlag = maximal lag between motor unit spike trains
%   jitter = tolerance in sec for the estimation of discharge times
%   tol = percentage of shared discharge times to define a duplicate
%   fsamp = sampling frequency

% Output:
%   Pulsenew = Pulse train of non-duplicated MU
%   distimenew = discharge times of non-duplicated MU
%   musclenew = indexes of non-duplicated MU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pulsenew, distimenew, musclenew, percom] = checkduplicates(PulseT, distime, muscle, maxlag, jitter, tol, fsamp)
f = waitbar(0,'Removing duplicates');
jit = round(jitter*fsamp);

% Generate binary spike trains
firings = zeros(size(PulseT));
for i = 1:size(PulseT,1)
    firings(i , distime{i}) = 1;
    distimmp{i} = [];
    for j = 1:jit
        distimmp{i} = [distimmp{i} distime{i}-j];
        distimmp{i} = [distimmp{i} distime{i}+j];
    end
    distimmp{i} = [distimmp{i} distime{i}];
end


MUn = length(distime);
x = 1/MUn;

i = 1;
% Remove duplicates
while ~isempty(distimmp)
    % Remove lag that may exist between MU
    distimetemp = cell(1,length(distimmp));
    for j = 1:length(distimmp)
        [c, lags] = xcorr(firings(1,:), firings(j,:), maxlag*2, 'coeff');
        [correl, idx] = max(c);
        if correl > 0.2
            distimetemp{j} = distimmp{j} + lags(idx);
        else
            distimetemp{j} = distimmp{j};
        end
    end
    
    % Find common discharge times
    comdis = zeros(1, length(distimmp));
    for j = 2:length(distimmp)
        com = intersect(distimmp{1}, distimetemp{j});
        com([false,diff(com) == 1]) = [];
        comdis(j) = length(com)/max([length(distime{1}) length(distime{j})]);
        clearvars com
    end
    
    % Flag duplicates and keep the MU with the lowest CoV of ISI
    duplicates = find(comdis >= tol);
    duplicates = [1 duplicates];
    percom{i} = comdis(duplicates);
    
    % Delete duplicates and save the surviving MU
    distimenew{i} = distime{duplicates(1)};
    Pulsenew(i,:) = PulseT(duplicates(1),:);
    musclenew(i,1) = muscle(duplicates(1));
    disp(duplicates)
    
    if length(duplicates) == 2
        musclenew(i,2) = muscle(duplicates(2));
    elseif length(duplicates) > 2
        tmpidx = setdiff(duplicates,duplicates(1));
        for j = 1:length(tmpidx)
            RT(j) = abs(distime{tmpidx(j)}(1)-distimenew{i}(1));
        end
        [~, survivor2] = min(RT);
        musclenew(i,2) = muscle(tmpidx(survivor2));
        tmpidx = setdiff(tmpidx, tmpidx(survivor2));
        duplicates = setdiff(duplicates,tmpidx);
        clearvars tmpidx RT
    else
        musclenew(i,2) = 0;
    end

    % Update firings and discharge times 
    for j = 1:length(duplicates)
        distimmp{duplicates(end-(j-1))} = [];
        distime{duplicates(end-(j-1))} = [];
    end
    distimmp = distimmp(~cellfun('isempty',distimmp));
    distime = distime(~cellfun('isempty',distime));

    firings(duplicates,:) = [];
    PulseT(duplicates,:) = [];
    muscle(duplicates) = [];
    waitbar(x*(MUn-length(distime)), f, [num2str(length(distime)) ' remaining MUs to check'])
    i = i + 1;
end

close(f);