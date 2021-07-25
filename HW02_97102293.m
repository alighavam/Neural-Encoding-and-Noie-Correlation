%% HW02 Advanced Neuroscience - Ali Ghavampour 97102293

%% Step 1 ===========================================================
clear all; close all; clc;
load UnitsData.mat;

% Plotting some of the units 
n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;

% psth = [];
% for i = 1:192
%     psth = [psth ; myPSTH(Unit,i,nbins,moveAVG)];
% end
t = linspace(-1.2,2,nbins);

% PSTH of  some units
unit_num = randi(n,[1,8]);
for i = 1:8
    num = unit_num(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(4,2,i)
    hold all
    plot(t,mean(psth,1),'k','LineWidth',2)
    xline(0,'--r','HandleVisibility','off');
    for cnd = 1:6
        plot(t,psth(cnd,:))
    end
    xlim([-1.2,2])
    legend("Avg","Cnd = 1","Cnd = 2","Cnd = 3","Cnd = 4","Cnd = 5","Cnd = 6")
    title("PSTH plot of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

%% Checking the vairability of units
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

freq_var = zeros(1,n);
hold all
for i = 1:n
    tmp = myPSTH(Unit,i,nbins,moveAVG);
    tmp = mean(tmp,1);
    freq_var(i) = max(tmp) - min(tmp);
    plot(t,tmp)
end
xline(-1,'--r','LineWidth',2);
xline(1.7,'--r','LineWidth',2);
title("Average PSTH over all conditions for each unit")
xlabel("Time(s)")
ylabel("Frequency(Hz)")
hold off

ind = find(freq_var>25);
figure;
% Choosing units by inspection
for i = 1:length(ind(1:10))
    num = ind(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(5,2,i)
    hold all
    plot(t,mean(psth,1),'k','LineWidth',2)
    xline(0,'--r','HandleVisibility','off');
    for cnd = 1:6
        plot(t,psth(cnd,:))
    end
    xlim([-1.2,2])
    legend("Avg","Cnd = 1","Cnd = 2","Cnd = 3","Cnd = 4","Cnd = 5","Cnd = 6")
    title("PSTH plot of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

%% Effect of reward and cue location on PSTH
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

unit_num = [12    30    35    55    69    73    91    96   119   134];
% unit_num = randi(n,[1,10]);
% based on reward
figure;
for i = 1:10
    num = unit_num(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(5,2,i)
    hold all
    plot(t,mean(psth(1:2,:),1),'b','LineWidth',2)
    plot(t,psth(1,:),'b','HandleVisibility','off')
    plot(t,psth(2,:),'b','HandleVisibility','off')

    plot(t,mean(psth(3:4,:),1),'r','LineWidth',2)
    plot(t,psth(3,:),'r','HandleVisibility','off')
    plot(t,psth(4,:),'r','HandleVisibility','off')

    plot(t,mean(psth(5:6,:),1),'k','LineWidth',2)
    plot(t,psth(5,:),'k','HandleVisibility','off')
    plot(t,psth(6,:),'k','HandleVisibility','off')
    
    legend("Reward = 3", "Reward = 6", "Reward = 9");
    xlim([-1.2,2]);
    title("PSTH of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end
% based on cue
figure;
for i = 1:10
    num = unit_num(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(5,2,i)
    hold all
    plot(t,mean(psth([1,3,5],:),1),'b','LineWidth',2)
    plot(t,psth(1,:),'b','HandleVisibility','off')
    plot(t,psth(3,:),'b','HandleVisibility','off')
    plot(t,psth(5,:),'b','HandleVisibility','off')

    plot(t,mean(psth([2,4,6],:),1),'r','LineWidth',2)
    plot(t,psth(2,:),'r','HandleVisibility','off')
    plot(t,psth(4,:),'r','HandleVisibility','off')
    plot(t,psth(6,:),'r','HandleVisibility','off')
    
    legend("Cue = -1", "Cue = 1");
    xlim([-1.2,2]);
    title("PSTH of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

%% Similar units psth plots
clear all; close all; clc;
load UnitsData.mat;

% Plotting some of the units 
n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 9;

t = linspace(-1.2,2,nbins);

% PSTH of  some units
unit_num = [3 12 52 66 70 79 111 217 255 274 279 418 419 440 463 480];
for i = 1:16
    num = unit_num(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(4,4,i)
    hold all
    plot(t,mean(psth,1),'k','LineWidth',2)
    xline(0,'--r','HandleVisibility','off');
    for cnd = 1:6
        plot(t,psth(cnd,:))
    end
    xlim([-1.2,2])
%     legend("Avg","Cnd = 1","Cnd = 2","Cnd = 3","Cnd = 4","Cnd = 5","Cnd = 6")
    title("PSTH plot of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

figure;
unit_num = [1 6 41 74 150 169 213 226 259 266 290 332 360 373 394 470];
for i = 1:16
    num = unit_num(i);
    [psth,~,Cond] = myPSTH(Unit,num,nbins,moveAVG);
    subplot(4,4,i)
    hold all
    plot(t,mean(psth,1),'k','LineWidth',2)
    xline(0,'--r','HandleVisibility','off');
    for cnd = 1:6
        plot(t,psth(cnd,:))
    end
    xlim([-1.2,2])
%     legend("Avg","Cnd = 1","Cnd = 2","Cnd = 3","Cnd = 4","Cnd = 5","Cnd = 6")
    title("PSTH plot of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

%% Acitivy analysis around t = 1.5ms
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60;
moveAVG = 0;
t = linspace(-1.2,2,nbins);

avgFiring = [];
avgFiringCND1 = [];
avgFiringCND2 = [];
avgFiringCND3 = [];
avgFiringCND4 = [];
avgFiringCND5 = [];
avgFiringCND6 = [];
for i = 1:n
    [~,psthAll,cnd] = myPSTH(Unit,i,nbins,moveAVG);
    for j = 1:192
        tmp = psthAll(j,:)/max(psthAll(j,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiring = [avgFiring, mean(tmp(48:53),2)];
        end
    end
    
    cnd1 = cnd{1};
    for j = 1:length(cnd1)
        ind = cnd1(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND1 = [avgFiringCND1, mean(tmp(48:53),2)];
        end
    end
    
    cnd2 = cnd{2};
    for j = 1:length(cnd2)
        ind = cnd2(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND2 = [avgFiringCND2, mean(tmp(48:53),2)];
        end
    end
    
    cnd3 = cnd{3};
    for j = 1:length(cnd1)
        ind = cnd3(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND3 = [avgFiringCND3, mean(tmp(48:53),2)];
        end
    end
    
    cnd4 = cnd{4};
    for j = 1:length(cnd4)
        ind = cnd4(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND4 = [avgFiringCND4, mean(tmp(48:53),2)];
        end
    end
    
    cnd5 = cnd{5};
    for j = 1:length(cnd5)
        ind = cnd5(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND5 = [avgFiringCND5, mean(tmp(48:53),2)];
        end
    end
    
    cnd6 = cnd{6};
    for j = 1:length(cnd1)
        ind = cnd6(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND6 = [avgFiringCND6, mean(tmp(48:53),2)];
        end
    end
end
bins = 50;
x = linspace(0,1,bins);
h = hist(avgFiring,x);
h1 = hist(avgFiringCND1,x);
h2 = hist(avgFiringCND2,x);
h3 = hist(avgFiringCND3,x);
h4 = hist(avgFiringCND4,x);
h5 = hist(avgFiringCND5,x);
h6 = hist(avgFiringCND6,x);

subplot(4,2,[1,2])
bar(x,h)
hold all
s = std(avgFiring);
m = mean(avgFiring);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")
hold off

subplot(4,2,3)
bar(x,h1)
hold on
s = std(avgFiringCND1);
m = mean(avgFiringCND1);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 1 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(4,2,4)
bar(x,h2)
hold on
s = std(avgFiringCND2);
m = mean(avgFiringCND2);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 2 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(4,2,5)
bar(x,h3)
hold on
s = std(avgFiringCND3);
m = mean(avgFiringCND3);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 3 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(4,2,6)
bar(x,h4)
hold on
s = std(avgFiringCND4);
m = mean(avgFiringCND4);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 4 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(4,2,7)
bar(x,h5)
hold on
s = std(avgFiringCND5);
m = mean(avgFiringCND5);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 5 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(4,2,8)
bar(x,h6)
hold on
s = std(avgFiringCND6);
m = mean(avgFiringCND6);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH condition 6 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

% Now for Rewards and Cue locations
avgFiring = [];
avgFiringRew3 = [];
avgFiringRew6 = [];
avgFiringRew9 = [];
avgFiringCue1 = []; % -1
avgFiringCue2 = []; % +1
for i = 1:n
    [~,psthAll,cnd] = myPSTH(Unit,i,nbins,moveAVG);
    for j = 1:192
        tmp = psthAll(j,:)/max(psthAll(j,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiring = [avgFiring, mean(tmp(48:53),2)];
        end
    end
    
    cnd1 = [cnd{1};cnd{2}];
    for j = 1:length(cnd1)
        ind = cnd1(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringRew3 = [avgFiringRew3, mean(tmp(48:53),2)];
        end
    end
    
    cnd2 = [cnd{3};cnd{4}];
    for j = 1:length(cnd2)
        ind = cnd2(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringRew6 = [avgFiringRew6, mean(tmp(48:53),2)];
        end
    end
    
    cnd3 = [cnd{5};cnd{6}];
    for j = 1:length(cnd1)
        ind = cnd3(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringRew9 = [avgFiringRew9, mean(tmp(48:53),2)];
        end
    end
    
    cnd4 = [cnd{1};cnd{3};cnd{5}];
    for j = 1:length(cnd4)
        ind = cnd4(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCue1 = [avgFiringCue1, mean(tmp(48:53),2)];
        end
    end
    
    cnd5 = [cnd{2};cnd{4};cnd{6}];
    for j = 1:length(cnd5)
        ind = cnd5(j);
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCue2 = [avgFiringCue2, mean(tmp(48:53),2)];
        end
    end
end
figure;
bins = 50;
x = linspace(0,1,bins);
h = hist(avgFiring,x);
h1 = hist(avgFiringRew3,x);
h2 = hist(avgFiringRew6,x);
h3 = hist(avgFiringRew9,x);
h4 = hist(avgFiringCue1,x);
h5 = hist(avgFiringCue2,x);

subplot(3,2,1)
bar(x,h)
hold on
s = std(avgFiring);
m = mean(avgFiring);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(3,2,2)
bar(x,h1)
hold on
s = std(avgFiringRew3);
m = mean(avgFiringRew3);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH Reward = 3 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(3,2,4)
bar(x,h2)
hold on
s = std(avgFiringRew6);
m = mean(avgFiringRew6);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH Reward = 6 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(3,2,6)
bar(x,h3)
hold on
s = std(avgFiringRew9);
m = mean(avgFiringRew9);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH Reward = 9 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(3,2,3)
bar(x,h4)
hold on
s = std(avgFiringCue1);
m = mean(avgFiringCue1);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH Cue = -1 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(3,2,5)
bar(x,h5)
hold on
s = std(avgFiringCue2);
m = mean(avgFiringCue2);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH Cue = +1 (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")


%% Activity analysis at different times
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60;
moveAVG = 9;
t = linspace(-1.2,2,nbins);

avgFiring = [];
avgFiringCND1 = [];
avgFiringCND2 = [];
avgFiringCND3 = [];
avgFiringCND4 = [];
avgFiringCND5 = [];
avgFiringCND6 = [];
for i = 1:n
    [~,psthAll,cnd] = myPSTH(Unit,i,nbins,moveAVG);
    for j = 1:192
        tmp = psthAll(j,:)/max(psthAll(j,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiring = [avgFiring, mean(tmp(20:25),2)];
        end
    end
    
    for j = 1:192
        ind = j;
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND1 = [avgFiringCND1, mean(tmp(48:53),2)];
        end
    end
    
    for j = 1:192
        ind = j;
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND2 = [avgFiringCND2, mean(tmp(1:6),2)];
        end
    end
    
    for j = 1:192
        ind = j;
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND3 = [avgFiringCND3, mean(tmp(55:60),2)];
        end
    end
    
    for j = 1:192
        ind = j;
        tmp = psthAll(ind,:)/max(psthAll(ind,:));
        nanfind = isnan(tmp);
        if isempty(find(nanfind))
            avgFiringCND4 = [avgFiringCND4, mean(tmp(38:43),2)];
        end
    end
end
bins = 50;
x = linspace(0,1,bins);
h = hist(avgFiring,x);
h1 = hist(avgFiringCND1,x);
h2 = hist(avgFiringCND2,x);
h3 = hist(avgFiringCND3,x);
h4 = hist(avgFiringCND4,x);

subplot(5,1,1)
bar(x,h2)
hold all
s = std(avgFiringCND2);
m = mean(avgFiringCND2);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around -1.2ms)")
ylabel("count")
xlabel("Normalized frequency")
hold off

subplot(5,1,2)
bar(x,h)
hold on
s = std(avgFiring);
m = mean(avgFiring);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 0ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(5,1,3)
bar(x,h4)
hold on
s = std(avgFiringCND4);
m = mean(avgFiringCND4);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 1ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(5,1,4)
bar(x,h1)
hold on
s = std(avgFiringCND1);
m = mean(avgFiringCND1);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 1.5ms)")
ylabel("count")
xlabel("Normalized frequency")

subplot(5,1,5)
bar(x,h3)
hold on
s = std(avgFiringCND3);
m = mean(avgFiringCND3);
xline(m-s/2,'k','LineWidth',2);
xline(m+s/2,'k','LineWidth',2);
xline(m,'r','LineWidth',2);
title("Histogram of PSTH frequency all conditions (around 2ms)")
ylabel("count")
xlabel("Normalized frequency")

%% Step 2 - GLM ============================================================
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60;
moveAVG = 0;
t = linspace(-1.2,2,nbins);

clc;
signif = [];
pValVec = [];
tic
for i = 1:n
    [~,psthAll,cnd] = myPSTH(Unit,i,nbins,moveAVG);
    y = zeros(length(cnd),1);
    reward = [3 3 6 6 9 9];
    for j = 1:6
        y(cnd{j}) = reward(j);
    end
    mdl = fitglm(psthAll,y);
    pVal = coefTest(mdl);
    if (pVal < 0.05)
        signif = [signif, i];
        pValVec = [pValVec, pVal];
    end
end
Test = [signif ; pValVec];
toc


%% Significant units PSTHs
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60;
moveAVG = 9;
t = linspace(-1.2,2,nbins);
unit_num = [12,13,15,17,50,52,53,72,109,163,165,172,176,247,272,291,297,347,...
    390,403,411,437,478];

for i = 1:16
    num = unit_num(i);
    psth = myPSTH(Unit,num,nbins,moveAVG);
    subplot(4,4,i)
    hold all
    plot(t,mean(psth,1),'k','LineWidth',2)
    xline(0,'--r','HandleVisibility','off');
    for cnd = 1:6
        plot(t,psth(cnd,:))
    end
    xlim([-1.2,2])
    title("PSTH plot of Unit " + num2str(num))
    xlabel("Time(s)")
    ylabel("Frequency(Hz)")
    hold off
end

%% Firing Rates Plots
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 9;
t = linspace(-1.2,2,nbins);

unit_num = [12,13,15,17,50,52,53,72,109,163,165,172,176,247,272,291,297,347,...
    390,403,411,437,478];

ind = randi(length(unit_num),[1,2])
NUM = unit_num(ind)

% NUM = [12,411]; % force the units
[~,psth1,cnd1] = myPSTH(Unit,NUM(1),nbins,moveAVG);
[~,psth2,cnd2] = myPSTH(Unit,NUM(2),nbins,moveAVG);

t_vec = [10,23,28,35,40,50,55,60];
for i = 1:8
    subplot(2,4,i)
    t = t_vec(i);
    x = psth1([cnd1{1}],t);
    y = psth2([cnd2{1}],t);
    scatter(x,y,'k','filled')
    hold on
    center = [mean(x),mean(y)];
    stdev = [std(x),std(y)];
    h = plotEllipses(center,stdev);
    h.EdgeColor = 'k'; 
    h.LineWidth = 2; 
    h.FaceColor = [0 0 0 .3]; %4th value is undocumented: transparency
    axis equal
    
    hold on
    x = psth1([cnd1{6}],t);
    y = psth2([cnd2{6}],t);
    scatter(x,y,'r','filled')
	legend("Condition 1", "Condition 6");
    hold on
    center = [mean(x),mean(y)];
    stdev = [std(x),std(y)];
    h = plotEllipses(center,stdev);
    h.EdgeColor = 'r'; 
    h.LineWidth = 2; 
    h.FaceColor = [0 0 1 .3]; %4th value is undocumented: transparency
    axis equal
    title("t = " + num2str(t*3.2/60-1.2))
    sgtitle("Units number "+num2str(NUM(1))+" and "+num2str(NUM(2)))
    xlabel("Unit "+num2str(NUM(1))+" frequency(Hz)")
    ylabel("Unit "+num2str(NUM(2))+" frequency(Hz)")
end

%% Step 3  ================================================================

%% PCA - Approach one
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

data = zeros(6,n,nbins);

for i = 1:n
    tmp = myPSTH(Unit,i,nbins,moveAVG);
    for cnd = 1:6
        data(cnd,i,:) = tmp(cnd,:);
    end
end

a1 = reshape(data(1,:,:),n,nbins)';
a1 = a1(22:end-1,:);
[eigvec1,projected1,eigval1,~,explained1,mu1] = pca(a1);
a2 = reshape(data(2,:,:),n,nbins)';
a2 = a2(22:end-1,:);
[eigvec2,projected2,eigval2,~,explained2,mu2] = pca(a2);
a3 = reshape(data(3,:,:),n,nbins)';
a3 = a3(22:end-1,:);
[eigvec3,projected3,eigval3,~,explained3,mu3] = pca(a3);
a4 = reshape(data(4,:,:),n,nbins)';
a4 = a4(22:end-1,:);
[eigvec4,projected4,eigval4,~,explained4,mu4] = pca(a4);
a5 = reshape(data(5,:,:),n,nbins)';
a5 = a5(22:end-1,:);
[eigvec5,projected5,eigval5,~,explained5,mu5] = pca(a5);
a6 = reshape(data(6,:,:),n,nbins)';
a6 = a6(22:end-1,:);
[eigvec6,projected6,eigval6,~,explained6,mu6] = pca(a6);

p = {projected1,projected2,projected3,projected4,projected5,projected6};
for i = 1:6
    subplot(3,2,i)
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot_dir(x,y);
    xlabel("PC 1");
    ylabel("PC 2");
    title("Condition " + num2str(i));
end

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 2");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");

figure
tmp = p{3};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2)
hold on
tmp = p{5};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{4};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{6};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2)

lgd = legend("Cue = left", "Cue = Right");
lgt.FontSize = 14;
xlabel("PC 1")
ylabel("PC 2")
zlabel("PC 3")
title("Visuallizing Different Cue Locations")

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,3);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 3");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");



%% Step 4 ================================================================

%% My shuffling approach
clear all; close all; clc;
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

dataShuffle = zeros(6,n,nbins);

for i = 1:n
    [~, psthAll,COND] = myPSTH(Unit,i,nbins,moveAVG);
    shuffled = psthAll(randperm(size(psthAll, 1)), :); % Shuffling rows
    for cnd = 1:6
        dataShuffle(cnd,i,:) = mean(shuffled(COND{cnd},:),1);
    end
end

% PCA
a1 = reshape(dataShuffle(1,:,:),n,nbins)';
a1 = a1(22:end-1,:);
[eigvec1,projected1,eigval1,~,explained1,mu1] = pca(a1);
a2 = reshape(dataShuffle(2,:,:),n,nbins)';
a2 = a2(22:end-1,:);
[eigvec2,projected2,eigval2,~,explained2,mu2] = pca(a2);
a3 = reshape(dataShuffle(3,:,:),n,nbins)';
a3 = a3(22:end-1,:);
[eigvec3,projected3,eigval3,~,explained3,mu3] = pca(a3);
a4 = reshape(dataShuffle(4,:,:),n,nbins)';
a4 = a4(22:end-1,:);
[eigvec4,projected4,eigval4,~,explained4,mu4] = pca(a4);
a5 = reshape(dataShuffle(5,:,:),n,nbins)';
a5 = a5(22:end-1,:);
[eigvec5,projected5,eigval5,~,explained5,mu5] = pca(a5);
a6 = reshape(dataShuffle(6,:,:),n,nbins)';
a6 = a6(22:end-1,:);
[eigvec6,projected6,eigval6,~,explained6,mu6] = pca(a6);

p = {projected1,projected2,projected3,projected4,projected5,projected6};
for i = 1:6
    subplot(3,2,i)
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot_dir(x,y);
    xlabel("PC 1");
    ylabel("PC 2");
    title("Condition " + num2str(i));
end

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 2");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");

figure
tmp = p{3};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2)
hold on
tmp = p{5};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{4};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{6};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2)

lgd = legend("Cue = left", "Cue = Right");
lgt.FontSize = 14;
xlabel("PC 1")
ylabel("PC 2")
zlabel("PC 3")
title("Visuallizing Different Cue Locations")

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,3);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 3");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");

%% CFR Shuffling
clear all; close all; clc;
startup
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

data = zeros(6,n,nbins);

for i = 1:n
    tmp = myPSTH(Unit,i,nbins,moveAVG);
    for cnd = 1:6
        data(cnd,i,:) = tmp(cnd,:);
    end
end
dataTensor = permute(data,[3 2 1]);

rng('shuffle','twister') % randomize the seed
surrogate_type = 'surrogate-TNC';
model_dim = 6;
times_msk = t>0 & t<2; % select movement-related times
[R2_data] = summarizeLDS(dataTensor(times_msk, :, :), model_dim, false); % function that evaluates the summary statistic of the LDS structure
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);
numSurrogates = 100;
params = [];
params.readout_mode = 2; % select readout mode (eg neuron mode)
params.shfl_mode = 3; % shuffle across tensor mode (eg condition mode)
params.fix_mode = 2; % shuffle per mode (shuffle for each neuron independently)
if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
else
    error('please specify a correct surrogate type') 
end
R2_surr = nan(numSurrogates, 1);
for i = 1:numSurrogates
    fprintf('surrogate %d from %d, ', i, numSurrogates)
    [surrTensor] = sampleCFR(dataTensor, params);       % generate CFR random surrogate data.
    [R2_surr(i)] = summarizeLDS(surrTensor(times_msk, :, :), model_dim, false);
end

P = mean(R2_data<= R2_surr); % (upper-tail test)

if P>=0.05
   fprintf('P value = %1.0e\n', P)
else
   fprintf('P value < %.3f\n', (P<0.001)*0.001 + (P<0.01 & P>=0.001)*0.01 + (P<0.05 & P>=0.01)*0.05)
end
%%%%%%%%%%%%%%%% plot null distribution
x = 0:0.03:1;
h = hist(R2_surr, x);
hf=figure;
set(hf, 'color', [1 1 1]);
hold on
box on
hb = bar(x, h);
set(hb,'facecolor',[0.5000    0.1059    0.0687],'barwidth',1,'edgecolor','none')
p = plot(R2_data, 0, 'ko', 'markerfacecolor', 'k', 'markersize',10);
xlabel('summary statistic (R^2)')
ylabel('count')
xlim([0 1])
set(gca, 'FontSize',12)
set(gca, 'xtick',[-1 0 1])
legend([p, hb], {'original data', surrogate_type})
legend boxoff

%% CFR Data PCA
clear all; close all; clc;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

load('CFRData.mat')
data = surrTensor;
data = permute(data,[3 2 1]);

% PCA
a1 = reshape(data(1,:,:),n,nbins)';
a1 = a1(22:end-1,:);
[eigvec1,projected1,eigval1,~,explained1,mu1] = pca(a1);
a2 = reshape(data(2,:,:),n,nbins)';
a2 = a2(22:end-1,:);
[eigvec2,projected2,eigval2,~,explained2,mu2] = pca(a2);
a3 = reshape(data(3,:,:),n,nbins)';
a3 = a3(22:end-1,:);
[eigvec3,projected3,eigval3,~,explained3,mu3] = pca(a3);
a4 = reshape(data(4,:,:),n,nbins)';
a4 = a4(22:end-1,:);
[eigvec4,projected4,eigval4,~,explained4,mu4] = pca(a4);
a5 = reshape(data(5,:,:),n,nbins)';
a5 = a5(22:end-1,:);
[eigvec5,projected5,eigval5,~,explained5,mu5] = pca(a5);
a6 = reshape(data(6,:,:),n,nbins)';
a6 = a6(22:end-1,:);
[eigvec6,projected6,eigval6,~,explained6,mu6] = pca(a6);

p = {projected1,projected2,projected3,projected4,projected5,projected6};
for i = 1:6
    subplot(3,2,i)
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot_dir(x,y);
    xlabel("PC 1");
    ylabel("PC 2");
    title("Condition " + num2str(i));
end

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 2");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");

figure
tmp = p{3};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2)
hold on
tmp = p{5};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{4};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{6};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2)

lgd = legend("Cue = left", "Cue = Right");
lgt.FontSize = 14;
xlabel("PC 1")
ylabel("PC 2")
zlabel("PC 3")
title("Visuallizing Different Cue Locations")

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,3);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 3");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");



%% TME Shuffling
clear all; close all; clc;
startup
load UnitsData.mat;

n = 481;
nbins = 60; % 3200 / 60 = 53.33ms bin sizes
moveAVG = 0;
t = linspace(-1.2,2,nbins);

data1 = zeros(6,n,nbins);

for i = 1:n
    tmp = myPSTH(Unit,i,nbins,moveAVG);
    for cnd = 1:6
        data1(cnd,i,:) = tmp(cnd,:);
    end
end
dataTensor = permute(data1,[3 2 1]);

rng('shuffle', 'twister') % randomize the seed
surrogate_type = 'surrogate-TNC';
model_dim = 6;
times_msk = t>0 & t<2; % select movement-related times
[R2_data] = summarizeLDS(dataTensor(times_msk, :, :), model_dim, false); % function that evaluates the summary statistic of the LDS structure
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);
numSurrogates = 100;
params = [];

if strcmp(surrogate_type, 'surrogate-T')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = [];
    params.margCov{3} = [];
    params.meanTensor = M.T;
elseif strcmp(surrogate_type, 'surrogate-TN')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = [];
    params.meanTensor = M.TN;
elseif strcmp(surrogate_type, 'surrogate-TNC')
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
else
    error('please specify a correct surrogate type') 
end

maxEntropy = fitMaxEntropy(params);             % fit the maximum entropy distribution
R2_surr = nan(numSurrogates, 1);
for i = 1:numSurrogates
    fprintf('surrogate %d from %d \n', i, numSurrogates)
    [surrTensor] = sampleTME(maxEntropy);       % generate TME random surrogate data.
    [R2_surr(i)] = summarizeLDS(surrTensor(times_msk, :, :), model_dim, false);
end

P = mean(R2_data<= R2_surr); % (upper-tail test)

if P>=0.05
   fprintf('P value = %1.0e\n', P)
else
   fprintf('P value < %.3f\n', (P<0.001)*0.001 + (P<0.01 & P>=0.001)*0.01 + (P<0.05 & P>=0.01)*0.05)
end
%%%%%%%%%%%%%%%% plot null distribution
x = 0:0.03:1;
h = hist(R2_surr, x);
hf=figure;
set(hf, 'color', [1 1 1]);
hold on
box on
hb = bar(x, h);
set(hb,'facecolor',[0.5000    0.3118    0.0176],'barwidth',1,'edgecolor','none')
p = plot(R2_data, 0, 'ko', 'markerfacecolor', 'k', 'markersize',10);
xlabel('summary statistic (R^2)')
ylabel('count')
xlim([0 1])
set(gca, 'FontSize',12)
set(gca, 'xtick',[-1 0 1])
legend([p, hb], {'original data', surrogate_type})
legend boxoff

data = surrTensor;
data = permute(data,[3 2 1]);

% PCA
a1 = reshape(data(1,:,:),n,nbins)';
a1 = a1(22:end-1,:);
[eigvec1,projected1,eigval1,~,explained1,mu1] = pca(a1);
a2 = reshape(data(2,:,:),n,nbins)';
a2 = a2(22:end-1,:);
[eigvec2,projected2,eigval2,~,explained2,mu2] = pca(a2);
a3 = reshape(data(3,:,:),n,nbins)';
a3 = a3(22:end-1,:);
[eigvec3,projected3,eigval3,~,explained3,mu3] = pca(a3);
a4 = reshape(data(4,:,:),n,nbins)';
a4 = a4(22:end-1,:);
[eigvec4,projected4,eigval4,~,explained4,mu4] = pca(a4);
a5 = reshape(data(5,:,:),n,nbins)';
a5 = a5(22:end-1,:);
[eigvec5,projected5,eigval5,~,explained5,mu5] = pca(a5);
a6 = reshape(data(6,:,:),n,nbins)';
a6 = a6(22:end-1,:);
[eigvec6,projected6,eigval6,~,explained6,mu6] = pca(a6);
figure;
p = {projected1,projected2,projected3,projected4,projected5,projected6};
for i = 1:6
    subplot(3,2,i)
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot_dir(x,y);
    xlabel("PC 1");
    ylabel("PC 2");
    title("Condition " + num2str(i));
end

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,2);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 2");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");

figure
tmp = p{3};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2)
hold on
tmp = p{5};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'k','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{4};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2,'HandleVisibility','off')
hold on
tmp = p{6};
x = tmp(:,1);
y = tmp(:,2);
z = tmp(:,3);
plot3(x,y,z,'r','LineWidth',2)

lgd = legend("Cue = left", "Cue = Right");
lgt.FontSize = 14;
xlabel("PC 1")
ylabel("PC 2")
zlabel("PC 3")
title("Visuallizing Different Cue Locations")

figure
for i = 1:6
    tmp = p{i};
    x = tmp(:,1);
    y = tmp(:,3);
    plot(x,y);
    hold on
    xlabel("PC 1");
    ylabel("PC 3");
    title("All Conditions");
end
legend("Cnd 1","Cnd 2","Cnd 3","Cnd 4","Cnd 5","Cnd 6");



%% Functions 
function [h1, h2] = plot_dir (vX, vY) % by Kangwon Lee
    %function [h1, h2] = plot_dir (vX, vY)
    %Plotting x-y variables with direction indicating vector to the next element.
    %Example
    %   vX = linspace(0,2*pi, 10)';
    %   vY = sin (vX);
    %   plot_dir(vX, vY);

    rMag = 0.5;

    % Length of vector
    lenTime = length(vX);

    % Indices of tails of arrows
    vSelect0 = 1:(lenTime-1);
    % Indices of tails of arrows
    vSelect1 = vSelect0 + 1;

    % X coordinates of tails of arrows
    vXQ0 = vX(vSelect0, 1);
    % Y coordinates of tails of arrows
    vYQ0 = vY(vSelect0, 1);

    % X coordinates of heads of arrows
    vXQ1 = vX(vSelect1, 1);
    % Y coordinates of heads of arrows
    vYQ1 = vY(vSelect1, 1);

    % vector difference between heads & tails
    vPx = (vXQ1 - vXQ0) * rMag;
    vPy = (vYQ1 - vYQ0) * rMag;

    % make plot 
    h1 = plot (vX, vY, '.-k','LineWidth',2); hold on;
    % add arrows 
    h2 = quiver (vXQ0,vYQ0, vPx, vPy, 0, 'r'); grid on; hold off
    set(h2,'AutoScale','on', 'AutoScaleFactor', 0.8)
    % axis equal
end

function [h1, h2] = plot_dir3 (vX, vY, vZ) % by Kangwon Lee
    %function [h1, h2] = plot_dir3 (vX, vY, vZ)
    %Plotting x-y variables with direction indicating vector to the next element.
    %Example
    %   vX = linspace(0,2*pi, 10)';
    %   vY = sin (vX);
    %   vZ = cos (vX);
    %   plot_dir3(vX, vY, vZ);

    rMag = 0.5;

    % Length of vector
    lenTime = length(vX);

    % Indices of tails of arrows
    vSelect0 = 1:(lenTime-1);
    % Indices of tails of arrows
    vSelect1 = vSelect0 + 1;

    % X coordinates of tails of arrows
    vXQ0 = vX(vSelect0, 1);
    % Y coordinates of tails of arrows
    vYQ0 = vY(vSelect0, 1);
    % X coordinates of tails of arrows
    vZQ0 = vZ(vSelect0, 1);

    % X coordinates of heads of arrows
    vXQ1 = vX(vSelect1, 1);
    % Y coordinates of heads of arrows
    vYQ1 = vY(vSelect1, 1);
    % Z coordinates of heads of arrows
    vZQ1 = vZ(vSelect1, 1);

    % vector difference between heads & tails
    vPx = (vXQ1 - vXQ0) * rMag;
    vPy = (vYQ1 - vYQ0) * rMag;
    vPz = (vZQ1 - vZQ0) * rMag;

    % make plot 
    h1 = plot3 (vX, vY, vZ, '.-k','LineWidth',2); hold on;
    % add arrows 
    h2 = quiver3 (vXQ0, vYQ0, vZQ0, vPx, vPy, vPz, 0, 'r'); grid on; hold off
    axis equal
end

function [P,vectors,values,M] = myPCA(A)
    M = mean(A,1);
    C = A - repmat(M,size(A,1),1);
    V = cov(C);
    values = eig(V);
    [vectors,~] = eig(V);
    P = vectors' * C';
    P = P';
end

function [x,y] = ellip(ellipse)
    x0 = ellipse.X0_in;
    y0 = ellipse.Y0_in;
    phi = ellipse.phi;
    a = ellipse.a;
    b = ellipse.b;
    t = 0 : 0.01 : 2*pi;
    x = x0 + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = y0 + b*sin(t)*cos(phi) + a*cos(t)*sin(phi);
end

function h = plotEllipses(cnt,rads,axh) % not written by me
    if nargin < 3 || isempty(axh)
       axh = gca();  
    end
    % Compute the lower, left corner of the rectangle.
    llc = cnt(:)-rads(:);
    % Compute the width and height
    wh = rads(:)*2; 
    % Draw rectangle 
    h = rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
end

function [psth,psthAll,cnd] = myPSTH(Unit,Unit_Num,nbins,moveAVG)
    trialsCell = Unit(Unit_Num).Trls;
%     nbins = 60; % 3200 / 60 = 53.33ms bins
    h = zeros(size(trialsCell,1),nbins);
    for i = 1:size(trialsCell,1)
        temp = trialsCell(i);
        temp = cell2mat(temp)';
        h(i,:) = hist(temp,nbins);
    end
    psthAll = h/(3.2/nbins);
    % Conditions
    Cnd =  Unit(12).Cnd;
    cnd1 = Cnd(1).TrialIdx;
    cnd2 = Cnd(2).TrialIdx;
    cnd3 = Cnd(3).TrialIdx;
    cnd4 = Cnd(4).TrialIdx;
    cnd5 = Cnd(5).TrialIdx;
    cnd6 = Cnd(6).TrialIdx;
    cnd = {cnd1,cnd2,cnd3,cnd4,cnd5,cnd6};

    % PSTH for 6 CNDs
    psth = zeros(6,nbins);
    for i = 1:6
        psth(i,:) = mean(h(cnd{i},:),1);
    end
    psth = psth/(3.2/nbins); % Changing unit to Hz
    if moveAVG ~= 0
        psth = movmean(psth,moveAVG,2);
        psthAll = movmean(psthAll,moveAVG,2);
    end
end
