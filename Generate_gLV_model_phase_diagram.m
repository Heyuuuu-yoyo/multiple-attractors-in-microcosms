%% Generate_gLV_model_phase_diagram
% gLV, S versus Alpha

close all;
clear;
T = 15000; % Simulation Iterative Time = T * step = 15000*0.1 = 1500
step = 0.1; % Simulation Step

Alpha = (60:-2:1)/20; % Species interaction strength (<\alpha_ij>)
SpeciesPool = 3:3:90; % Species pool size (S)
ReplicateNumber = 256; % Communities number every pixel

num_physical_cores = feature('numcores');
if isempty(gcp('nocreate'))
    parpool('local', 64, 'IdleTimeout', 120);
end

FlucMatrix = zeros(length(Alpha),length(SpeciesPool));
MultiMatrix = zeros(length(Alpha),length(SpeciesPool));
StatesMatrix = zeros(length(Alpha),length(SpeciesPool));
GlobalMatrix = zeros(length(Alpha),length(SpeciesPool));

FlucCell = cell(length(Alpha),length(SpeciesPool)); % Fluctuation states number
StableCell = cell(length(Alpha),length(SpeciesPool)); % Total Stable states number
TotalCell = cell(length(Alpha),length(SpeciesPool)); 

DSSCell = cell(length(Alpha),length(SpeciesPool)); % Different Stable states Number
FlucSSCell = cell(length(Alpha),length(SpeciesPool)); % Different Fluc states Number
GlobalCell = cell(length(Alpha),length(SpeciesPool));

StatesLabelCell = cell(length(Alpha),length(SpeciesPool));
StableDiversityCell = cell(length(Alpha),length(SpeciesPool)); 
MeanStableDiversityCell = cell(length(Alpha),length(SpeciesPool)); 

tempStatesLabel = cell(ReplicateNumber, 1);
tempStatesAbundance = cell(ReplicateNumber, 1);
tempStableDiversity = cell(ReplicateNumber, 1);

% Multistable threshold in silico
threshold = 0.01;

for k = 1:length(Alpha)
    meanA = Alpha(k); % gLV Interaction Strength
    for j = 1:length(SpeciesPool)
        S = SpeciesPool(j); % Species number
        Time = 100; % Different initial condition number
        r = ones(1,S); % Intrinsic Growth Rate
        InitCon = ones(S,S) * 1e-4 + eye(S,S)*0.1;% S different initial conditions(species abundances)
        StableStatesnumber = [];
        MeanDiversity = [];
        Fluc = [];
        Global = [];
        TotalStableStates = [];
        parfor rep = 1:ReplicateNumber
            AA = 2*meanA*rand(S,S); % Uniform Distribution
%             AA = (meanA/sqrt(3))*randn(S,S) + meanA; % Gaussian Distribution
            for i =1:S
                AA(i,i) = 1;
            end

            localStateLabels = [];
            localStableDiversity = [];

            Std_total = [];
            States = [];
            for time = 1:Time
                N = zeros(T,S);
                if (time<=S)
                    N(1,:) = InitCon(time,:);
                else
                    N(1,:) = rand(1,S)/S*2;
                end
                for i = 2:T
                    k1 = step* N(i-1,:).*(r' - AA*N(i-1,:)')';
                    k2 = step* (N(i-1,:)+k1/2).*(r' - AA*(N(i-1,:)+k1/2)')';
                    k3 = step* (N(i-1,:)+k2/2).*(r' - AA*(N(i-1,:)+k2/2)')';
                    k4 = step* (N(i-1,:)+k3).*(r' - AA*(N(i-1,:)+k3)')';
                    N(i,:) = N(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6 + 10^-6*ones(1,S)*step;
                end

                if ((max(max(abs(N))))<5)  % To ensure the community has converged
                    Std = std(N(T-1000:T,:));
                    States = [States;mean(N(T-1000:T,:))];
                    Std_total = [Std_total,max(Std)];
                end
            end
            if (size(States,1)>0)
                Stable_states = States(Std_total<0.001,:);
                Stable_states_relative = (Stable_states'./sum(Stable_states'))';
            else
                Stable_states = [];
                Stable_states_relative = [];
            end
            if (size(Stable_states,1)>0)
                num_rows = size(Stable_states_relative, 1);
                state_labels = zeros(num_rows, 1); 
                current_state = 0; 
                states_temp = []; 
                for i = 1:num_rows
                    if state_labels(i) == 0 
                        current_state = current_state + 1;
                        state_labels(i) = current_state;
                        states_temp = [states_temp; Stable_states_relative(i, :)]; 
                        for i2 = i+1:num_rows
                            if state_labels(i2) == 0
                                if all(abs(Stable_states_relative(i2, :) - Stable_states_relative(i, :)) < threshold)
                                    state_labels(i2) = current_state; 
                                end
                            end
                        end
                    end
                end
                localStateLabels = state_labels;
                localStableDiversity = sum(Stable_states'>1e-3)';
                num_states = current_state;
                StableStatesnumber = [StableStatesnumber, num_states];
                MeanDiversity = [MeanDiversity, mean(sum(Stable_states'>1e-3))]; % Diversity
            else
                StableStatesnumber = [StableStatesnumber, 0];
                MeanDiversity = [MeanDiversity, 0]; % Diversity
                num_states = 0;
            end
            TotalStableStates = [TotalStableStates,size(Stable_states,1)];
            if(size(Stable_states,1)<size(States,1))
                Fluc = [Fluc, size(States,1)-size(Stable_states,1)];
            else
                Fluc = [Fluc, 0];
            end
            if(size(Stable_states,1)==size(States,1) && num_states==1)
                Global = [Global,1];
            else
                Global = [Global,0];
            end
            tempStatesLabel{rep} = localStateLabels;
            tempStableDiversity{rep} = localStableDiversity;
            N = [];
        end

        StatesMatrix(k,j) = mean(StableStatesnumber);
        MultiMatrix(k,j) = mean(StableStatesnumber>1);
        GlobalMatrix(k,j) = mean(Global);
        FlucMatrix(k,j) = mean(Fluc>0);
        
        FlucCell{k,j} = Fluc;
        StableCell{k,j} = TotalStableStates;
        TotalCell{k,j} = Fluc + TotalStableStates;
        
        DSSCell{k,j} = StableStatesnumber;
        GlobalCell{k,j} = Global;
        MeanStableDiversityCell{k,j} = MeanDiversity;

        StatesLabelCell{k, j} = tempStatesLabel;
        StableDiversityCell{k, j} = tempStableDiversity;
    end
end

a = length(Alpha);
b = length(SpeciesPool);
GlobalMatrix = zeros(a, b); % Global stability
OnlyFlucMatrix = zeros(a, b); % Pure Fluctuation 
OnlyMultiStableMatrix = zeros(a, b); % Pure multistable
MultiAttractorMatrix = zeros(a, b); % Multiple attractors
MultiStable_Matrix = zeros(a, b); % Has multiple stable states
FlucStable_Matrix = zeros(a, b); % Fluctuation-Stable
Fluc_Matrix = zeros(a, b); % Has fluctuating states

TOTAL = ReplicateNumber;
for i = 1:a
    for j = 1:b
        GlobalMatrix(i, j) = sum(GlobalCell{i,j}==1)/TOTAL;
        OnlyFlucMatrix(i, j) = sum(StableCell{i,j}==0 & FlucCell{i,j}>0)/TOTAL;
        MultiAttractorMatrix(i, j) = sum(DSSCell{i,j}>0 & FlucCell{i,j}>0)/TOTAL + sum(DSSCell{i,j}>1 & FlucCell{i,j}==0)/TOTAL;
        OnlyMultiStableMatrix(i,j) = sum(DSSCell{i,j}>1 & FlucCell{i,j}==0)/TOTAL;
        MultiStable_Matrix(i,j) = sum(DSSCell{i,j}>1)/TOTAL;
        FlucStable_Matrix(i,j) = sum(DSSCell{i,j}>0 & FlucCell{i,j}>0)/TOTAL;
        Fluc_Matrix(i, j) = sum(FlucCell{i,j}>0)/TOTAL;
    end
end
