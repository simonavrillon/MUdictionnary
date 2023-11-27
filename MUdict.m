%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUDict.m
% To run this script, add to the path the functions bandpassingals.m,
% getMUfilters.m, getPulseT.m, checkduplicate.m, and extend.m
%
% Then, paste in the current folder the files from the same participant
% (e.g., from S1_10_DF.otb+_decomp.mat_edited.mat to
% S1_80_DF.otb+_decomp.mat_edited.mat)
%
% Output:
% RTm = Recruitment threshold average across intensities
% RT = Recruitment threshold across intensities
% Disrate = Time stamps, instantaneous discharge rates, force, and target
% for all the intensities where the motor unit was tracked
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: Get the filters from intensity 1 and apply them on itensity 2

files = dir('*edited.mat');
for n = 1:size(files,1)-1
    data1 = load(files(n).name); % Load file 1
    data2 = load(files(n+1).name); % Load file 2
    nmu1 = 1;
    nmu2 = 1;
    PulseT1 = [];
    PulseT2 = [];
    Distime1 = {};
    Distime2 = {};
    
    for j = 1:size(data1.edition.Pulsetrain,2) % Loop on each grid
        if size(data1.edition.Pulsetrain{j},1) > 0
            EMGmask = data1.signal.EMGmask{j} + data2.signal.EMGmask{j}; % Get the masks for boths grids (we remove all the channels masked on either task 1 or task 2
            for k = 1:size(data1.edition.Pulsetrain{j},1) % Format the discharge times for the 1st grid
                Distim{k} = data1.edition.Dischargetimes{j,k};
                idxmu1(nmu1) = j*100+k; % Use an index for each motor unit (e.g., 101 is the first motor unit of grid 1)
                idx = find(data1.signal.target > 0);
                Distime = intersect(idx, data1.edition.Dischargetimes{j,k});
                DRt = 1./(diff((Distime/data1.signal.fsamp)));
                for l = 1:length(DRt)
                    timeDR(l) = round((Distime(l+1)-Distime(l))/2 + Distime(l)) / data1.signal.fsamp;
                    to(l) = data1.signal.path(round((Distime(l+1)-Distime(l))/2 + Distime(l)));
                    ta(l) = data1.signal.target(round((Distime(l+1)-Distime(l))/2 + Distime(l)));
                end
                
                Disrate{n,nmu1} = cell(1,8);
                Disrate{n,nmu1}{n} = [timeDR', DRt', to', ta'];
        
                tmp = find(DRt > 5 & DRt < 25, 1 ,'first');
                if ~isempty(tmp)
                    firingprop{n}(nmu1,1) = mean(data1.signal.path(round(timeDR(tmp)*data1.signal.fsamp)-12:round(timeDR(tmp)*data1.signal.fsamp)+12));
                end
        
                firingprop{n}(nmu1,2) = idxmu1(nmu1);
    
                RTf{n,nmu1} = zeros(2,8);
                RTf{n,nmu1}(1,n) = firingprop{n}(nmu1,1);
                RTf{n,nmu1}(2,n) = 1;
    
                clearvars Distime timeDR DRt idx ta to
                nmu1 = nmu1 + 1;
            end
            if ~isempty(k)
                datatmp1 = data1.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp1 = bandpassingals(datatmp1,data1.signal.fsamp, 1);
                datatmp2 = data2.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp2 = bandpassingals(datatmp2,data2.signal.fsamp, 1);
    
                MUFilters = getMUfilters(datatmp1, EMGmask, Distim); % Get the MU filters from the first file
                [PulseT, Distime, ~] = getPulseT(datatmp2, EMGmask, MUFilters, data1.signal.fsamp);  % Apply the MU filters on the second file
                PulseT1 = [PulseT1; PulseT]; % Concatenate all the grids (file 1)
                Distime1 = [Distime1, Distime];
            end
            clearvars Distime PulseT Distim MUFilters EMGmask datatmp1 datatmp2
        end
    end

    for j = 1:size(data2.edition.Pulsetrain,2) % Loop on each grid
        if size(data2.edition.Pulsetrain{j}, 1) > 0
            PulseT2 = [PulseT2; data2.edition.Pulsetrain{j}]; % Concatenate all the grids (file 2)
            for k = 1:size(data2.edition.Pulsetrain{j},1)
                Distime2{nmu2} = data2.edition.Dischargetimes{j,k};
                idxmu2(nmu2) = j*1000+k; % Use an index for each motor unit (e.g., 1001 is the first motor unit of grid 1)
                nmu2 = nmu2+1;
            end
        end
    end

    Pulsetemp = [PulseT2; PulseT1]; % Concatenate the Pulse trains generated with the filters of the file 1 and the Pulse trains from file 2
    Distemp = [Distime2 Distime1];
    idxmus = [idxmu2 idxmu1];
    [Pulsetemp, Distemp, matchetemp] = checkduplicates(Pulsetemp, Distemp, idxmus, round(data2.signal.fsamp/40), 0.00025, 0.3, data2.signal.fsamp); % Find the matches
    
    % Here you will need to sort the matches and the new Pulse trains
    % 2 values in a raw (matches between files 1 and 2)
    % 1 value in the hundreds (New pulse train without a match), if editable,
    % new MU and new match
    % 1 value in the thousands, no match
    [~, idx] = sort(matchetemp(:,2), 'descend');
    
    for i = 1:length(idx)
        matches{n}(i,:) = matchetemp(idx(i), :);
        PulseT{n}(i,:) = Pulsetemp(idx(i), :);
        distimenew{n}{i} = Distemp{idx(i)};
    end
    clearvars -except files n matches firingprop RTf Disrate
end

nmu1 = 1;
n = size(files,1);
data1 = load(files(n).name); % Load file 1
for j = 1:size(data1.edition.Pulsetrain,2) % Loop on eeach grid
    if size(data1.edition.Pulsetrain{j},1) > 0
        EMGmask = data1.signal.EMGmask{j}; % Get the masks for boths grids (we remove all the channels masked on either task 1 or task 2
        for k = 1:size(data1.edition.Pulsetrain{j},1) % Format the discharge times for the 1st grid
            Distim{k} = data1.edition.Dischargetimes{j,k};
            idxmu1(nmu1) = j*100+k; % Use an index for each motor unit (e.g., 101 is the first motor unit of grid 1)
            idx = find(data1.signal.target > 0);
            Distime = intersect(idx, data1.edition.Dischargetimes{j,k});
            DRt = 1./(diff((Distime/data1.signal.fsamp)));
            for l = 1:length(DRt)
                timeDR(l) = round((Distime(l+1)-Distime(l))/2 + Distime(l)) / data1.signal.fsamp;
                to(l) = data1.signal.path(round((Distime(l+1)-Distime(l))/2 + Distime(l)));
                ta(l) = data1.signal.target(round((Distime(l+1)-Distime(l))/2 + Distime(l)));
            end
            
            Disrate{n,nmu1} = cell(1,8);
            Disrate{n,nmu1}{n} = [timeDR', DRt', to', ta'];
            
            tmp = find(DRt > 5 & DRt < 25, 1 ,'first');
            if ~isempty(tmp)
                firingprop{n}(nmu1,1) = mean(data1.signal.path(round(timeDR(tmp)*data1.signal.fsamp)-12:round(timeDR(tmp)*data1.signal.fsamp)+12));
            end
    
            firingprop{n}(nmu1,2) = idxmu1(nmu1);
    
            RTf{n,nmu1} = zeros(2,8);
            RTf{n,nmu1}(1,n) = firingprop{n}(nmu1,1);
            RTf{n,nmu1}(2,n) = 1;
        
            clearvars Distime timeDR DRt idx ta to
            nmu1 = nmu1 + 1;
        end
    end
end
clearvars -except files n matches firingprop RTf Disrate

%% Step 2: Get the filters from intensity 1 and apply them on itensity 3

files = dir('*edited.mat');
for n = 1:size(files,1)-2
    data1 = load(files(n).name); % Load file 1
    data2 = load(files(n+2).name); % Load file 2
    nmu1 = 1;
    nmu2 = 1;
    PulseT1 = [];
    PulseT2 = [];
    Distime1 = {};
    Distime2 = {};
    for j = 1:size(data1.edition.Pulsetrain,2) % Loop on each grid
        if size(data1.edition.Pulsetrain{j},1) > 0
            EMGmask = data1.signal.EMGmask{j} + data2.signal.EMGmask{j}; % Get the masks for boths grids (we remove all the channels masked on either task 1 or task 2
            for k = 1:size(data1.edition.Pulsetrain{j},1) % Format the discharge times for the 1st grid
                Distim{k} = data1.edition.Dischargetimes{j,k};
                idxmu1(nmu1) = j*100+k; % Use an index for each motor unit (e.g., 101 is the first motor unit of grid 1)
                nmu1 = nmu1 + 1;
            end
            if ~isempty(k)
                datatmp1 = data1.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp1 = bandpassingals(datatmp1,data1.signal.fsamp, 1);
                datatmp2 = data2.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp2 = bandpassingals(datatmp2,data2.signal.fsamp, 1);
    
                MUFilters = getMUfilters(datatmp1, EMGmask, Distim); % Get the MU filters from the first file
                [PulseT, Distime, ~] = getPulseT(datatmp2, EMGmask, MUFilters, data1.signal.fsamp);  % Apply the MU filters on the second file
                PulseT1 = [PulseT1; PulseT]; % Concatenate all the grids (file 1)
                Distime1 = [Distime1, Distime];
            end
            clearvars Distime PulseT Distim MUFilters EMGmask datatmp1 datatmp2
        end
    end

    for j = 1:size(data2.edition.Pulsetrain,2) % Loop on each grid
        if size(data2.edition.Pulsetrain{j},1) > 0
            PulseT2 = [PulseT2; data2.edition.Pulsetrain{j}]; % Concatenate all the grids (file 2)
            for k = 1:size(data2.edition.Pulsetrain{j},1)
                Distime2{nmu2} = data2.edition.Dischargetimes{j,k};
                idxmu2(nmu2) = j*1000+k; % Use an index for each motor unit (e.g., 1001 is the first motor unit of grid 1)
                nmu2 = nmu2+1;
            end
        end
    end

    Pulsetemp = [PulseT2; PulseT1]; % Concatenate the Pulse trains generated with the filters of the file 1 and the Pulse trains from file 2
    Distemp = [Distime2 Distime1];
    idxmus = [idxmu2 idxmu1];
    [Pulsetemp, Distemp, matchetemp] = checkduplicates(Pulsetemp, Distemp, idxmus, round(data2.signal.fsamp/40), 0.00025, 0.3, data2.signal.fsamp); % Find the matches
    
    % Here you will need to sort the matches and the new Pulse trains
    % 2 values in a raw (matches between files 1 and 2)
    % 1 value in the hundreds (New pulse train without a match), if editable,
    % new MU and new match
    % 1 value in the thousands, no match
    [~, idx] = sort(matchetemp(:,2), 'descend');
    
    for i = 1:length(idx)
        matches{2,n}(i,:) = matchetemp(idx(i), :);
        PulseT{n}(i,:) = Pulsetemp(idx(i), :);
        distimenew{n}{i} = Distemp{idx(i)};
    end
    clearvars -except files n matches firingprop RTf Disrate
end

%% Step 3: Reorganise the data according to the tracking
% Here we have three catagories: The unmatches with i+1 and i+2, the unmatches with i+1 
% but match with i+2, and the matches with i+1
% 1 identify the unmatches

for n = 1:size(files,1)-2
    % Identify the unmatches with i+1 but match with i+2
    [~,IA,IB] = intersect(matches{1,n}(:,1), matches{2,n}(:,2));

    if ~isempty(IB)
        matches{1,n}(IA,:) = [];
        % Rearrange the idxs
        matches{2,n}(IB,3) = floor(matches{2,n}(IB,1)/1000);
        matches{2,n}(IB,4) = matches{2,n}(IB,1) - matches{2,n}(IB,3)*1000 + matches{2,n}(IB,3)*100;
        matches{2,n}(matches{2,n}(:,4) == 0,:) = [];
        
        % Identify the values for the matched units
        [~,IA,IB] = intersect(matches{2,n}(:,2), firingprop{n}(:,2));
        for i = 1:length(IB)
            temp(i) = matches{2,n}(IA(i),4);
            tempRTf{i} = RTf{n,IB(i)};
            for j = 1:n
                tempDisrate{i,j} = Disrate{n,IB(i)}{j};
            end
        end
        
        [~,IA,IB] = intersect(temp, firingprop{n+2}(:,2));
        % Update for the next iteration
        for i = 1:length(IB)
            RTf{n+2,IB(i)} = RTf{n+2,IB(i)} + tempRTf{IA(i)};
            for j = 1:n
                Disrate{n+2,IB(i)}{j} = tempDisrate{IA(i),j};
            end
        end
        clearvars temp tempRTf tempDisrate
    end
end

nmu = 1;
for n = 1:size(files,1)-1
    % Identify the unmatches with both i+1 and i+2
    [~,~,IB] = intersect(matches{1,n}(:,1), firingprop{n}(:,2));
    for i = 1:length(IB)
        RTfc(nmu,:) = RTf{n,IB(i)}(1,:);
        Disratec{nmu} = Disrate{n,IB(i)};
        nmu = nmu+1;
    end
    clearvars IB

    % Rearrange the idxs
    matches{1,n}(:,3) = floor(matches{1,n}(:,1)/1000);
    matches{1,n}(:,4) = matches{1,n}(:,1) - matches{1,n}(:,3)*1000 + matches{1,n}(:,3)*100;

    % Identify the matches with i+1
    [~,IA,IB] = intersect(matches{1,n}(:,2), firingprop{n}(:,2));
    if ~isempty(IB)
        for i = 1:length(IB)
            temp(i) = matches{1,n}(IA(i),4);
            tempRTf{i} = RTf{n,IB(i)};
            for j = 1:n
                tempDisrate{i,j} = Disrate{n,IB(i)}{j};
            end
        end
        
        [C,IA,IB] = intersect(temp, firingprop{n+1}(:,2));
        % Update for the next iteration
        for i = 1:length(IB)
            RTf{n+1,IB(i)} = RTf{n+1,IB(i)} + tempRTf{IA(i)};
            for j = 1:n
                Disrate{n+1,IB(i)}{j} = tempDisrate{IA(i),j};
            end
        end
        clearvars temp tempRTf tempDisrate
    end
end

for i = 1:size(firingprop{size(files,1)},1)
    RTfc(nmu,:) = RTf{size(files,1),i}(1,:);
    Disratec{nmu} = Disrate{size(files,1),i};
    nmu = nmu+1;
end

RTfc(RTfc==0) = NaN;

%% Step 4: Inspect the RTs to split the MUs wrongly matched (with >10%MVC between RTs)

iter = size(RTfc,1);
flag = ones(1,iter);
while sum(flag) > 0
iter = size(RTfc,1);
flag = zeros(1,iter);
for i = 1:iter
    tmp = diff(RTfc(i,:));
    if max(abs(tmp)) > 10
        flag(i) = 1;
        idx = find(abs(tmp) > 10, 1, 'first');
        RTfc(end+1,1:8) = NaN;
        RTfc(end,1:idx) = RTfc(i,1:idx);
        RTfc(end+1,1:8) = NaN;
        RTfc(end,idx+1:end) = RTfc(i,idx+1:end);
        
        Disratec{end+1} = cell(1,8);
        for j = 1:idx
            Disratec{end}{j} = Disratec{i}{j};
        end
        Disratec{end+1} = cell(1,8);
        for j = idx+1:8
            Disratec{end}{j} = Disratec{i}{j};
        end
        flag(end+2) = 0;
    else
        flag(i) = 0;
    end
end

RTfc(flag==1,:) = [];

nmu = 1;
for i = 1:length(flag)
    if flag (i) == 0
        Disratecc{nmu} = Disratec{i};
        nmu = nmu + 1;
    end
end
clearvars MUAPc Disratec
Disratec = Disratecc;
clearvars MUAPcc Disratecc
end

clearvars -except RTfc Disratec
%% Step 5: Clean the data and save the output
RTtemp = nanmean(RTfc,2);
[~,idx] = sort(RTtemp,'ascend');

for i = 1:length(idx)
    RT(i,:) = RTfc(idx(i),:);
    Disrate{i} = Disratec{idx(i)};
end
RTm = nanmean(RT,2);

clearvars RTtemp RTfc idx i Disratec