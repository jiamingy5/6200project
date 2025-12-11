function mef_sensitivity_simple
clc; clear; close all;

%% ========== 1. Read Data ==========
flowP  = readtable('inflow_P.csv','VariableNamingRule','preserve');
unregP = readtable('unreginflow_P.csv','VariableNamingRule','preserve');
evapP  = readtable('evap_P.csv','VariableNamingRule','preserve');
elevP  = readtable('elevation_P.csv','VariableNamingRule','preserve');
elevM  = readtable('elevation_M.csv','VariableNamingRule','preserve');

flowP.datetime  = datetime(flowP.datetime);
unregP.datetime = datetime(unregP.datetime);
evapP.datetime  = datetime(evapP.datetime);

data = outerjoin(flowP, unregP, "Keys","datetime","MergeKeys",true);
data = outerjoin(data, evapP, "Keys","datetime","MergeKeys",true);
data = sortrows(data,"datetime");

numVars = varfun(@isnumeric, data, "OutputFormat", "uniform");
data{:,numVars} = fillmissing(data{:,numVars},"constant",0);

vnames = data.Properties.VariableNames;
idxUnreg  = find(contains(vnames,"unreg"),1);
idxEvap   = find(contains(vnames,"evap"),1);

%% ========== 2. Inflow(unregulated inflow) ==========
Q_cfs = data{:,idxUnreg};
bad = (Q_cfs < 0) | (Q_cfs > 8e4);
Q_cfs(bad) = NaN;
Q_cfs = fillmissing(Q_cfs,'linear','EndValues','nearest');

Q = Q_cfs * 86400 / 43560;   % acre-ft/day
EvapP = data{:,idxEvap};

fprintf('\n===== Cleaned Q stats (acre-ft/day) =====\n');
fprintf('min = %.0f\n', min(Q));
fprintf('mean = %.0f\n', mean(Q));
fprintf('max = %.0f\n\n', max(Q));

%% ========== 3. Stage–storage ==========
Powell_SS = [
    3525  5.93e6
    3575  9.52e6
    3666  1.929e7
    3700  2.432e7];

P_elev = Powell_SS(:,1); 
P_stor = Powell_SS(:,2);

PowellElev = @(S) interp1(P_stor, P_elev, S,'pchip','extrap');
PowellStor = @(H) interp1(P_elev, P_stor, H,'pchip','extrap');

Mead_SS = [
    950   2.006e6
    1050  7.683e6
    1205  2.3936e7
    1220  2.6399e7
    1229  2.762e7];

M_elev = Mead_SS(:,1);
M_stor = Mead_SS(:,2);

MeadElev = @(S) interp1(M_stor, M_elev, S,'pchip','extrap');
MeadStor = @(H) interp1(M_elev, M_stor, H,'pchip','extrap');

%% ========== 4. Initial Storage Based on Starting Elevations ==========
colHP = elevP.Properties.VariableNames{contains(elevP.Properties.VariableNames,"pool") | contains(elevP.Properties.VariableNames,"elev")};
colHM = elevM.Properties.VariableNames{contains(elevM.Properties.VariableNames,"pool") | contains(elevM.Properties.VariableNames,"elev")};

SP0 = PowellStor(elevP.(colHP)(1));
SM0 = MeadStor(elevM.(colHM)(1));

%% ========== 5. Model Parameters ==========
minElevP = 3525;                  
Rbase    = 20548;    % Baseline release (acre-ft/day)
DM       = 7500000/365;   % Annual Lower Basin delivery requirement → daily

%% ========== 6. Realistic MEF Range for Sensitivity Testing ==========
% 0–10,000 acre-ft/day ≈ 0–5000 cfs
MEFs = (0:2000:10000)';    
nM = numel(MEFs);

p_sat    = zeros(nM,1);
p_lowP   = zeros(nM,1);
tierFrac = zeros(nM,3);

tier1 = 1075;
tier2 = 1050;

fprintf('=== Running simplified MEF sensitivity ===\n');

for i = 1:nM
    MEF = MEFs(i);

    [R,HP,HM] = simulate(MEF, SP0, SM0, Q, EvapP, ...
                         PowellStor, PowellElev, ...
                         MeadStor, MeadElev, ...
                         DM, Rbase);

    p_sat(i)  = mean(R >= MEF);
    p_lowP(i) = mean(HP <= minElevP);

    tierFrac(i,1) = mean(HM >= tier1);
    tierFrac(i,2) = mean(HM < tier1 & HM >= tier2);
    tierFrac(i,3) = mean(HM < tier2);
end

%% ========== 7. Plot Results ==========
figure;
plot(MEFs, p_sat,'o-','LineWidth',2);
grid on;
title('MEF Satisfaction Probability');
xlabel('MEF (acre-ft/day)');
ylabel('P(R ≥ MEF)');

figure;
plot(MEFs, p_lowP,'s-','LineWidth',2);
grid on;
title('Powell Low-Elevation Risk');
xlabel('MEF (acre-ft/day)');
ylabel('P(Elev ≤ 3525 ft)');

figure;
plot(MEFs, tierFrac(:,1),'^-','LineWidth',2); hold on;
plot(MEFs, tierFrac(:,2),'v-','LineWidth',2);
plot(MEFs, tierFrac(:,3),'d-','LineWidth',2);
grid on;
title('Mead Tier Distribution vs MEF');
xlabel('MEF (acre-ft/day)');
ylabel('Fraction of Time');
legend('Tier 1','Tier 2','Tier 3','Location','best');

fprintf('=== MEF sensitivity finished ===\n');

end


%% ========== Simplified Powell–Mead Simulation Model ==========
function [R, HP, HM] = simulate(MEF, SP0, SM0, Q, EvapP, ...
                                PowellStor, PowellElev, ...
                                MeadStor, MeadElev, ...
                                DM, Rbase)

T = numel(Q);
SP = SP0;
SM = SM0;

R  = zeros(T,1);
HP = zeros(T,1);
HM = zeros(T,1);

for t = 1:T

    Rtarget = max(Rbase, MEF);      
    avail   = max(SP + Q(t) - EvapP(t), 0);
    R(t)    = min(avail, Rtarget);

    SP = SP + Q(t) - R(t) - EvapP(t);
    SP = max(SP,0);

    SM = SM + R(t) - DM;
    SM = max(SM,0);

    HP(t) = PowellElev(SP);
    HM(t) = MeadElev(SM);
end
end
