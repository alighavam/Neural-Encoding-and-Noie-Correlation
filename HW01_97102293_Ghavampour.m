%% HW01 Advanced Neuroscience - Ali Ghavampour 97102293

%% Question01 =========================================================
clear all; close all; clc;

% Part A 
fr = 100;
tSim = 1;
nTrials = 1000;
[spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);

% Part B
count = sum(spikes,2);
x = 60:140;
h = hist(count,x,'k')/nTrials;
subplot(2,2,1)
bar(x,h)
hold on
plot(x, poisspdf(x,fr),'r','LineWidth',1.5);
title("Spike Count Probability")
xlabel("lambda")
ylabel("Probability")


% Part C
binSize = 1;
x = [1:binSize:100];
h = zeros(nTrials,size(x,2));
PoissonTau = [];
for i = 1:nTrials
    spikeTime = find(spikes(i,:));
    spikeIntervals = spikeTime(2:length(spikeTime)) - spikeTime(1:length(spikeTime) - 1);
    PoissonTau = [PoissonTau , spikeIntervals]; % Saving all the ISIs
    h(i,:) = hist(spikeIntervals, x);
end
CV_PoissonSpikes = std(PoissonTau)/mean(PoissonTau) % CV of Spikes
h = sum(h,1)/100/nTrials;
subplot(2,2,2)
bar(x, h);
t = 0:100;
y = fr * exp(-fr * t/1000)/1000;
hold on
plot(t,y,'r','LineWidth',1.5);
title("ISIs Histogram Alongside Theoretical")
xlabel("tau(ms)")
ylabel("Probability")

% Repeating Part A to C while removing every K spike
k = 4;
newSpikes = renewalProcess(spikes, k);
count = sum(newSpikes,2);
x = 15:35;
h = hist(count,x,'k')/nTrials;
subplot(2,2,3)
bar(x,h)
title("Spike Count Probability After Removing Procedure")
xlabel("lambda")
ylabel("Probability")

binSize = 1;
x = [1:binSize:100];
h = zeros(nTrials,size(x,2));
NewTau = [];
for i = 1:nTrials
    spikeTime = find(newSpikes(i,:));
    spikeIntervals = spikeTime(2:length(spikeTime)) - spikeTime(1:length(spikeTime) - 1);
    NewTau = [NewTau , spikeIntervals]; % Saving All the ISIs
    h(i,:) = hist(spikeIntervals, x);
end
CV_DeletedSpikes = std(NewTau)/mean(NewTau) % CV of Spikes After Procedure
h = sum(h,1)/(fr/4)/nTrials;
subplot(2,2,4)
bar(x,h);
title("ISIs Histogram After Removing Procedure")
xlabel("tau(ms)")
ylabel("Probability")

%% Part G - Refractory Period
clear all; close all; clc;

% Theoretical Reconstruction
figure;
delta_t = 1:0.1:30;
CV1 = 1/sqrt(1) * ((delta_t - 1)./delta_t);
CV4 = 1/sqrt(4) * ((delta_t - 1)./delta_t);
CV51 = 1/sqrt(51) * ((delta_t - 1)./delta_t);
scatter(delta_t, CV1, 'k')
hold on
scatter(delta_t, CV4, 'k')
hold on
scatter(delta_t, CV51, 'k')
hold on
yline(0.5,'k');
hold on
yline(1,'k');
hold on
yline(0.14,'k');
ylim([0 1.1])
xlabel('$$\overline{\Delta t}$$','interpreter','latex');
ylabel("CV")
title("Figure 6 of Paper, Theoretical Reconstruction")

% Reconstruction based on simulation
tic
f = [999:-50:200,199:-10:150,150:-2:30];
CV1 = [];
for fr = f
    tSim = 10;
    nTrials = 500;
    [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
    binSize = 1;
    x = [1:binSize:100];
    h = zeros(nTrials,size(x,2));
    PoissonTau = [];
    for i = 1:nTrials
        spikeTime = find(spikes(i,:));
        spikeIntervals = spikeTime(2:length(spikeTime)) - spikeTime(1:length(spikeTime) - 1);
        PoissonTau = [PoissonTau , spikeIntervals]; % Saving all the ISIs
        h(i,:) = hist(spikeIntervals, x);
    end
    PoissonTau = PoissonTau + 1; % The refractory period is 1ms
    CV_PoissonSpikes = std(PoissonTau)/mean(PoissonTau); % CV of Spikes
    CV1 = [CV1, CV_PoissonSpikes];
end
figure;
scatter(1000./f, CV1, 'k')
hold on

CV4 = [];
for fr = f
    tSim = 10;
    nTrials = 500;
    [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
    spikes = renewalProcess(spikes, 4);
    binSize = 1;
    x = [1:binSize:100];
    h = zeros(nTrials,size(x,2));
    PoissonTau = [];
    for i = 1:nTrials
        spikeTime = find(spikes(i,:));
        spikeIntervals = spikeTime(2:length(spikeTime)) - spikeTime(1:length(spikeTime) - 1);
        PoissonTau = [PoissonTau , spikeIntervals]; % Saving all the ISIs
        h(i,:) = hist(spikeIntervals, x);
    end
    PoissonTau = PoissonTau + 1; % The refractory period is 1ms
    CV_PoissonSpikes = std(PoissonTau)/mean(PoissonTau); % CV of Spikes
    CV4 = [CV4, CV_PoissonSpikes];
end
scatter(1000./f, CV4, 'k')
hold on

CV51 = [];
for fr = f
    tSim = 10;
    nTrials = 500;
    [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
    spikes = renewalProcess(spikes, 51);
    binSize = 1;
    x = [1:binSize:100];
    h = zeros(nTrials,size(x,2));
    PoissonTau = [];
    for i = 1:nTrials
        spikeTime = find(spikes(i,:));
        spikeIntervals = spikeTime(2:length(spikeTime)) - spikeTime(1:length(spikeTime) - 1);
        PoissonTau = [PoissonTau , spikeIntervals]; % Saving all the ISIs
        h(i,:) = hist(spikeIntervals, x);
    end
    PoissonTau = PoissonTau + 1; % The refractory period is 1ms
    CV_PoissonSpikes = std(PoissonTau)/mean(PoissonTau); % CV of Spikes
    CV51 = [CV51, CV_PoissonSpikes];
end
scatter(1000./f, CV51, 'k')
hold on
yline(0.5,'k');
hold on
yline(1,'k');
hold on
yline(0.14,'k');
ylim([0 1.1])
xlabel('$$\overline{\Delta t}$$','interpreter','latex');
ylabel("CV")
title("Figure 6 of Paper, Simulation Reconstruction")
toc

%% Question02 ==========================================================

%% Part A
clear all; close all; clc;

RI = 20*10^-3;
taum = 13*10^-3;
vth = 15*10^-3;
vr = 0;
tsim = 100*10^-3;
dt = 0.001;
t = 0:dt:tsim-dt;
v = zeros(size(t));
v(1) = vr;
for i = 2:length(t)
    f = (-v(i-1) + RI)/taum;
    v(i) = v(i-1) + dt*f; 
    if (v(i) > vth)
        v(i) = 0;
    end
end
plot(t*1000,v*1000,'k')
xlabel("t(ms)")
ylabel("V(mv)")
title("Simulation of LIF model for 100ms")

%% Part C - Building and Testing The Model
clear all; close all; clc;

% Parameters
dt = 1/10000;
R = 0.001;
vth = 15*10^-3;
vr = 0;
taum = 13*10^-3; % 13ms
tr = 0.001; % refractory period
tau_peak = 1.5*10^-3; % 1.5ms

% Spike Train
fr = 200;
tSim = 1;
nTrials = 1;
[spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
InputRate = length(find(spikes))

% Kernel
t1 = 0:dt:15/1000-dt;
k = t1 .* exp(-t1/(tau_peak));
k = k / max(k);

% Making I(t)
I_exc = 20*conv(spikes,k,'same');
subplot(2,1,1)
plot(tVec,I_exc,'k')
hold on
plot(tVec,spikes*10,'r')
title("Neual Input I(t)")
xlabel("t(ms)");
ylabel("I(t)")
lgd = legend("I(t)","Input Spikes");
lgd.FontSize = 12;

[v, spikes] = LIF(vr,R,tr,taum,dt,vth,I_exc,tVec);
OutputRate = length(find(spikes))
x = find(spikes);
isi = (x(2:end)-x(1:end-1))*dt*1000;
subplot(2,1,2)
v(find(spikes==5)) = 50*10^-3;
plot(tVec*1000,v*1000,'k');
hold on;
plot(tVec*1000,spikes,'r');
title("Simulation of LIF Neuron")
xlabel("t(ms)");
ylabel("v(mV)")
lgd = legend("V(t)","Output Spikes");
lgd.FontSize = 12;
        
%% Part C - Reconstructing Figure 8
clear all; close all; clc;

% Parameters
dt = 1/10000;
R = 0.001;
vth = 15*10^-3;
vr = 0;
taum = 13*10^-3; % 13ms
tr = 0.001; % refractory period
tau_peak = 1.5*10^-3; % 1.5ms

% Spike Train
fr = 750;
tSim = 1;
nTrials = 1;
[spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);

taum_mat = [0.1:0.1:1 , 2:11] * 10^-3;
k_mat = [1:10,20:10:100];
f = [72,94,100,105,131,138,143,149,156,163,230,290,340,384,430,474,520,570,610,660];

CV = zeros(length(k_mat),length(taum_mat));
tic
for j = 1:length(k_mat)
    k = k_mat(j);
    for i = 1:length(taum_mat)
        taum = taum_mat(i);
        fr = f(i);
        isi = [];
        for iter = 1:15
            tSim = 10;
            nTrials = 1;
            [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
            t1 = 0:dt:15/1000-dt;
            kernel = t1 .* exp(-t1/(tau_peak));
            kernel = kernel / max(kernel);
            I_exc = 20*conv(spikes,kernel,'same');
            [v, spikes] = LIF(vr,R,tr,taum,dt,vth,I_exc,tVec);
            spikes = renewalProcess(spikes, k);
            x = find(spikes);
            isi = [isi,(x(2:end)-x(1:end-1))*dt*1000];
        end
        CV(j,i) = std(isi)/mean(isi);
    end
end
toc

[row,col] = find(0.05<CV & CV<0.15);
loglog(taum_mat(col)*1000,k_mat(row),'k');
hold on
[row,col] = find(0.15<CV & CV<0.25);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.25<CV & CV<0.35);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.35<CV & CV<0.45);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.45<CV & CV<0.55);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.55<CV & CV<0.65);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.65<CV & CV<0.75);
loglog(taum_mat(col)*1000,k_mat(row));
hold on
[row,col] = find(0.75<CV & CV<0.85);
loglog(taum_mat(col)*1000,k_mat(row));
hold on

lgd = legend("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8");
lgd.FontSize = 10;

%% Part C - Effect of width and magnitude of EPSCs on Cv
clear all; close all; clc;

% Parameters
dt = 1/10000;
R = 0.001;
vth = 15*10^-3;
vr = 0;
taum = 13*10^-3; % 13ms
tr = 0.001; % refractory period
tau_peak = 1.5*10^-3; % 1.5ms
tSim = 1;
fr = 200;

mag_mat = [1:5:2000];
tpeak_mat = [0.1:0.1:50] * 10^-3;

cv = zeros(size(mag_mat));
tic
for i = 1:length(mag_mat)
    mag = mag_mat(i);
    isi = [];
    for j = 1:15
        [spikes, tVec] = poissonSpikeGen(fr, tSim, 1, 0);
        I = mag * inputCur(dt,tau_peak,spikes);
        [v, spikes] = LIF(vr,R,tr,taum,dt,vth,I,tVec);
        temp = find(spikes);
        isi = [isi, (temp(2:end) - temp(1:end-1))*dt*1000];
    end
    cv(i) = std(isi)/mean(isi);
end
toc
subplot(2,1,1)
scatter(mag_mat,cv,'k');
title("Effect of EPSCs Magnitude on CV");
xlabel("Coeficient of EPSCs")
ylabel("CV")

cv = zeros(size(tpeak_mat));
tic
for i = 1:length(tpeak_mat)
    tau_peak = tpeak_mat(i);
    isi = [];
    for j = 1:15
        [spikes, tVec] = poissonSpikeGen(fr, tSim, 1, 0);
        I = 100 * inputCur(dt,tau_peak,spikes);
        [v, spikes] = LIF(vr,R,tr,taum,dt,vth,I,tVec);
        temp = find(spikes);
        isi = [isi, (temp(2:end) - temp(1:end-1))*dt*1000];
    end
    cv(i) = std(isi)/mean(isi);
end
toc
subplot(2,1,2)
scatter(tpeak_mat,cv,'k');
ylim([0.4 1]);
title("Effect of EPSCs Width on CV");
xlabel("Tau Peak")
ylabel("CV")




%% Part D - Inhibitory neurons
clear all; close all; clc;

% Parameters
dt = 1/10000;
R = 0.001;
vth = 15*10^-3;
vr = 0;
taum = 13*10^-3; % 13ms
tr = 0.001; % refractory period
tau_peak = 1.5*10^-3; % 1.5ms

% making I(t) mat
neurons = 500;
percents = [50:2:100]/100;
cv = zeros(size(percents));
tic
for j = 1:length(percents)
    exc = percents(j);
    inh = 1-exc;
    I_exc = [];
    tic
    for i = 1:neurons*exc
        fr = 100;
        tSim = 1;
        nTrials = 1;
        [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 0);
        I_exc = [I_exc ; inputCur(dt,tau_peak,spikes)];
    end
    I_inh = [];
    for i = 1:neurons*inh
        fr = 100;
        tSim = 1;
        nTrials = 1;
        [spikes, tVec] = poissonSpikeGen(fr, tSim, nTrials, 1);
        I_inh = [I_inh ; inputCur(dt,tau_peak,spikes)];
    end
    if (inh == 0)
        I_inh = 0;
    end
    toc
    I = sum(I_exc,1) + sum(I_inh,1);
    [v, spikes] = LIF(vr,R,tr,taum,dt,vth,I,tVec);
    temp = find(spikes);
    isi = (temp(2:end) - temp(1:end-1))*dt*1000;
    std(isi)/mean(isi)
    cv(j) = std(isi)/mean(isi);
end
toc
plot((1-percents)*100,cv,'k');
title("CV of a Network of 500 Excitatory and Inhibitory Neurons");
xlabel("Percentage of Inhibitory Neurons");
ylabel("CV");

%% Part E 
clear all; close all; clc;
M = 200;
ratio_mat = 0.1:0.01:1;
D_mat = [5:1:100]*10^-3;
dt = 1/100000;

[spikes, tVec] = poissonSpikeGen(10, 1, M, 0);
cv = zeros(length(ratio_mat),length(D_mat));

tic
for i = 1:length(ratio_mat)
    ratio = ratio_mat(i);
    N = M * ratio;
    output = zeros(size(tVec));
    for j = 1:length(D_mat)
        D = D_mat(j);
        ind = round(D / dt);
        for k = 1:ind:length(tVec)-ind+1
            tmp = spikes(:,k:k+ind-1);
            tmp = sum(tmp,2);
            x = length(find(tmp));
            if (x >= N)
                output(k+ind-1) = 1;
            end
        end
        tmp = find(output);
        isi = (tmp(2:end)-tmp(1:end-1))*dt*1000;
        cv(i,j) = std(isi)/mean(isi);
        std(isi)/mean(isi);
    end
end
toc

CV = cv;
% [row,col] = find(0.08<CV & CV<0.12);
% loglog(D_mat(col)*1000,ratio_mat(row),'k');
% hold on
[row,col] = find(0.18<CV & CV<0.22);
scatter(D_mat(col)*1000,ratio_mat(row),[],[151,151,151]./255,'filled');
hold on
[row,col] = find(0.28<CV & CV<0.32);
scatter(D_mat(col)*1000,ratio_mat(row),'k','filled');
hold on
[row,col] = find(0.38<CV & CV<0.42);
scatter(D_mat(col)*1000,ratio_mat(row),[],[46,222,245]./255,'filled');
hold on
[row,col] = find(0.48<CV & CV<0.52);
scatter(D_mat(col)*1000,ratio_mat(row),'g','filled');
hold on
[row,col] = find(0.58<CV & CV<0.62);
scatter(D_mat(col)*1000,ratio_mat(row),'r','filled');
hold on
[row,col] = find(0.68<CV & CV<0.72);
scatter(D_mat(col)*1000,ratio_mat(row),'m','filled');
hold on
[row,col] = find(0.78<CV & CV<0.82);
scatter(D_mat(col)*1000,ratio_mat(row),'b','filled');

set(gca, 'XScale', 'log')
title("Effect of N/M and D(window) on CV");
xlabel("D(ms)");
ylabel("N/M");
lgd = legend("0.2","0.3","0.4","0.5","0.6","0.7","0.8");
lgd.FontSize = 14;

%% Part F 
clear all; close all; clc;

neurons = 200;
exc = 0.8*neurons;
inh = 0.2*neurons;
threshold_mat = 1:exc;
D_mat = [5:1:100]*10^-3;
dt = 1/100000;

[spikesExc, tVec] = poissonSpikeGen(10, 1, exc, 0);
[spikesInh, tVec] = poissonSpikeGen(10, 1, inh, 0);
cv = zeros(length(threshold_mat),length(D_mat));

tic
for i = 1:length(threshold_mat)
    th = threshold_mat(i);
    output = zeros(size(tVec));
    for j = 1:length(D_mat)
        D = D_mat(j);
        ind = round(D / dt);
        for k = 1:ind:length(tVec)-ind+1
            tmpExc = spikesExc(:,k:k+ind-1);
            tmpExc = sum(tmpExc,2);
            x1 = length(find(tmpExc));
            
            tmpInh = spikesInh(:,k:k+ind-1);
            tmpInh = sum(tmpInh,2);
            x2 = length(find(tmpInh));
            
            if (x1-x2 >= th)
                output(k+ind-1) = 1;
            end
        end
        tmp = find(output);
        isi = (tmp(2:end)-tmp(1:end-1))*dt*1000;
        cv(i,j) = std(isi)/mean(isi);
        std(isi)/mean(isi);
    end
end
toc

CV = cv;
% [row,col] = find(0.08<CV & CV<0.12);
% loglog(D_mat(col)*1000,ratio_mat(row),'k');
% hold on
[row,col] = find(0.18<CV & CV<0.22);
scatter(D_mat(col)*1000,threshold_mat(row),[],[151,151,151]./255,'filled');
hold on
[row,col] = find(0.28<CV & CV<0.32);
scatter(D_mat(col)*1000,threshold_mat(row),'k','filled');
hold on
[row,col] = find(0.38<CV & CV<0.42);
scatter(D_mat(col)*1000,threshold_mat(row),[],[46,222,245]./255,'filled');
hold on
[row,col] = find(0.48<CV & CV<0.52);
scatter(D_mat(col)*1000,threshold_mat(row),'g','filled');
hold on
[row,col] = find(0.58<CV & CV<0.62);
scatter(D_mat(col)*1000,threshold_mat(row),'r','filled');
hold on
[row,col] = find(0.68<CV & CV<0.72);
scatter(D_mat(col)*1000,threshold_mat(row),'m','filled');
hold on
[row,col] = find(0.78<CV & CV<0.82);
scatter(D_mat(col)*1000,threshold_mat(row),'b','filled');

set(gca, 'XScale', 'log')
title("Effect of threshold and D(window) on CV");
xlabel("D(ms)");
ylabel("Threshold");
lgd = legend("0.2","0.3","0.4","0.5","0.6","0.7","0.8");
lgd.FontSize = 14;


%% Functions 
function [v, spikes] = LIF(vr,R,tr,taum,dt,vth,I,tVec)
    v = zeros(size(tVec));
    spikes = zeros(size(tVec));
    v(1) = vr;
    i = 2;
    while (i <= length(tVec))
        temp = v(i-1) + dt * (-v(i-1) + R*I(i-1))/taum;
        if (temp < vth)
            v(i) = temp;
            i = i+1;
        else
            spikes(i) = 5;
            ind = i + tr/dt;
            v(i:ind) = vr;
            i = ind+1;
        end
    end
    v(find(spikes==5)) = 50*10^-3;
end

function I = inputCur(dt,tau_peak,spikeTrain)
    spikes = spikeTrain;
    t1 = 0:dt:15/1000-dt;
    k = t1 .* exp(-t1/(tau_peak));
    I = conv(spikes,k,'same');
    I = I/max(abs(I));
end

function renewdProcess = renewalProcess(spikes, k) % removing every kth spike
    renewdProcess = zeros(size(spikes));
    for i = 1:size(spikes,1)
        temp = find(spikes(i,:));
        indexes = temp(1:k:length(temp));
        renewdProcess(i,indexes) = 1;
    end
end

function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials, inhibit) 
    dt = 1/10000; % s
    nBins = floor(tSim/dt); 
    spikeMat = rand(nTrials, nBins) < fr*dt; 
    if inhibit == 1
        spikeMat = -1 * spikeMat;
    end
    tVec = 0:dt:tSim-dt;
end

