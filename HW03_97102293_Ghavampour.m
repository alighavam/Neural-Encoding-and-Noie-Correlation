%% Advance Neuro HW03 - Ali Ghavampour - 97102293
clear all; close all; clc;

%% Questions 01 =====================================================
clear all; close all; clc;

binWidth = 20;
nMonkey = 3;
t = 0:binWidth:1000-binWidth;
gratsNum = 12;
monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
IND = {'IND_monkey1.mat','IND_monkey2.mat','IND_monkey3.mat'};
bestNeurons = zeros(nMonkey,gratsNum,3);

for imonkey = 1:nMonkey
    load(monkeys{imonkey});
    load(IND{imonkey});
    keepNeurons = out;
    n = size(S(1).mean_FRs,1);
    ind = find(keepNeurons);
    figure;
    for igrat = 1:gratsNum
        psth = S(igrat).mean_FRs;
        meanVec = mean(psth,2);
        indMax1 = find(meanVec == max(meanVec)); % First Maximum
        indMax2 = find(meanVec == max(meanVec(meanVec<max(meanVec)))); % Second Maximum
        thirdMax = max(meanVec((meanVec<max(meanVec(meanVec<max(meanVec))))));
        indMax3 = find(meanVec == thirdMax);
        bestNeurons(imonkey,igrat,1) = ind(indMax1);
        bestNeurons(imonkey,igrat,2) = ind(indMax2);
        bestNeurons(imonkey,igrat,3) = ind(indMax3);
        subplot(6,2,igrat)
        for i = 1:n
            plot(t,psth(i,:),'k','HandleVisibility','off')
            if (i == indMax1)
                h1 = plot(t,psth(i,:),'r','LineWidth',1.8,'DisplayName',sprintf('Neuron %d',ind(i)));
                yline(max(meanVec),'--r','HandleVisibility','off');
            end
            if (i == indMax2)
                h2 = plot(t,psth(i,:),'b','LineWidth',1.8,'DisplayName',sprintf('Neuron %d',ind(i))');
                yline(max(meanVec(meanVec<max(meanVec))),'--b','HandleVisibility','off');
            end
%             if (i == indMax3)
%                 h3 = plot(t,psth(i,:),'g','LineWidth',1.8,'DisplayName',sprintf('Neuron %d',ind(i))');
%                 yline(thirdMax,'--g','HandleVisibility','off');
%             end
            hold on
        end
        legend([h1,h2]);
        title("PSTH of Neurons - Grating " + num2str(igrat))
        xlabel("t(ms)")
        ylabel("firing rate")
    end
    sgtitle(sprintf("Monkey %d",imonkey));
end

%% Raster Plot
clear all; close all; clc;

monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
nMonkey = 3;
neuron = 25;
grat = [1,2,3,4,5,6];
load(monkeys{1});
for j = 1:6
    trials = zeros(200,1000);
    for i = 1:200
        trials(i,:) = S(grat(j)).trial(i).spikes(neuron,:);
    end
    subplot(3,2,j)
    RasterPlot(trials, 1000)
    title(sprintf("grating %d",grat(j)))
    xlabel("Time(ms)")
    ylabel("Trial")
end
sgtitle(sprintf("Monkey 1 - Neuron %d",neuron))

%% Tuning Curves
clear all; close all; clc;
binWidth = 20;
nMonkey = 3;
t = 0:binWidth:1000-binWidth;
gratsNum = 12;
monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
IND = {'IND_monkey1.mat','IND_monkey2.mat','IND_monkey3.mat'};
gratX = 0:30:330;

neurons = [35, 45, 97]; % Chosen neurons for each monkey
tuningCurve = zeros(nMonkey,gratsNum);
for imonkey = 1:nMonkey
    load(monkeys{imonkey});
    load(IND{imonkey});
    keepNeurons = out;
    activeNeurons = find(keepNeurons);
    ind = find(activeNeurons == neurons(imonkey));
    for igrat = 1:gratsNum
        tmp = S(igrat).mean_FRs;
        tuningCurve(imonkey,igrat) = mean(tmp(ind,:));
    end
end
for i = 1:3
    subplot(1,3,i)
    plot(gratX, tuningCurve(i,:), 'k', 'LineWidth', 2)
    ylabel("mean activity")
    xlabel("preferred stimulus")
    title(sprintf("Tuning Curve of Neuron %d of Monkey %d",neurons(i),i));
    xlim([0,330])
    grid on
end


%% Question 02 ===========================================================
clear all; close all; clc;
binWidth = 20;
nMonkey = 3;
t = 0:binWidth:1000-binWidth;
gratsNum = 12;
monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
IND = {'IND_monkey1.mat','IND_monkey2.mat','IND_monkey3.mat'};
dataS = {'data_monkey1_gratings.mat','data_monkey2_gratings.mat','data_monkey3_gratings.mat'};

colorMat = zeros(3,10,10);
for imonkey = 1:nMonkey
    load(dataS{imonkey});
    load(IND{imonkey});
    load(monkeys{imonkey});
    map = data.MAP;
    removeInd = find(out == 0);
    ch = data.CHANNELS;
    ch(removeInd,:) = [];
    ch(:,2) = [];
    map(isnan(map)) = 0;
    % for i = 1:length(removeInd)
    %     tmp = removeInd(i);
    %     map(find(map == tmp)) = 0;
    % end
    for row = 1:size(map,1)
        for col = 1:size(map,2)
            unit = map(row,col);
            if (unit ~= 0 && length(find(ch == unit)) ~= 0)
                neuron = find(ch == unit);
                neuron = neuron(1);
                tuningCurve = zeros(1,6);
                for igrat = 1:6
                    tmp = S(igrat).mean_FRs;
                    tmp2 = S(igrat+6).mean_FRs;
                    tuningCurve(igrat) = (mean(tmp(neuron,:)) + mean(tmp2(neuron,:)))/2;
%                     tuningCurve(igrat) = mean(tmp(neuron,:));
                end
                prefered = find(tuningCurve == max(tuningCurve));
%                 if (prefered > 6)
%                     prefered = prefered - 6;
%                 end
                colorMat(imonkey,row,col) = prefered;
            end
        end
    end
end

x = 0:0.4:3.6;
for i = 1:3
    subplot(1,3,i)
    colormap(jet)
%     s = pcolor(flip(reshape(colorMat(i,:,:),10,10)));
    imagesc(x,x,reshape(colorMat(i,:,:),10,10));
%     s.FaceColor = 'interp';
    title(sprintf("Monkey %d",i))
end

figure
for i = 1:3
    subplot(1,3,i)
    colormap(jet)
    pcolor(x,x,flip(reshape(colorMat(i,:,:),10,10)));
    shading interp;
    title(sprintf("Monkey %d",i))
end

%% Question 03 - Figure A ==========================================================
clear all; close all; clc;
binWidth = 20;
nMonkey = 3;
t = 0:binWidth:1000-binWidth;
gratsNum = 12;
monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
IND = {'IND_monkey1.mat','IND_monkey2.mat','IND_monkey3.mat'};
dataS = {'data_monkey1_gratings.mat','data_monkey2_gratings.mat','data_monkey3_gratings.mat'};

r1 = {};
r2 = {};
r3 = {};
r4 = {};
for imonkey = 1:nMonkey
    rtmp1 = [0;0];
    rtmp2 = [0;0];
    rtmp3 = [0;0];
    rtmp4 = [0;0];
    tic
    load(dataS{imonkey});
    load(IND{imonkey});
    load(monkeys{imonkey});
    toc
    map = data.MAP;
    removeInd = find(out == 0);
    ch = data.CHANNELS;
    ch(removeInd,:) = [];
    ch(:,2) = [];
    map(isnan(map)) = 0;

    tunMat = zeros(12,length(find(out)));
    for i = 1:length(find(out))
        tuningCurve = zeros(1,12);
        for igrat = 1:12
            tmp = S(igrat).mean_FRs;
            tuningCurve(igrat) = mean(tmp(i,:));
        end
        tunMat(:,i) = tuningCurve;
    end
    corrMat = corr(tunMat);

    tic
    for i = 1:length(find(out))-1
        for j = i+1:length(find(out))
            rSig = corrMat(i,j);
            if (rSig > 0.5)
                % Finding r_sc
                tunCur01 = tunMat(:,i);
                tunCur02 = tunMat(:,j);
                pref01 = find(tunCur01 == max(tunCur01));
                pref02 = find(tunCur02 == max(tunCur02));
                count01 = zeros(200,1);
                count02 = zeros(200,2);
                for trl = 1:200
                    tmp01 = S(pref01).trial(trl).counts;
                    tmp02 = S(pref02).trial(trl).counts;
                    count01(trl) = sum(tmp01(i,:));
                    count02(trl) = sum(tmp02(j,:));
                end
                % ====================================================!!!!!!!!!!!!!!!!!
                zsc = zscore([count01, count02]);
                rsc = corr(zsc);
                rsc = rsc(2);

                % finding distance
                unit01 = ch(i);
                unit02 = ch(j);
                if (unit01 ~= unit02)
                    [row01, col01] = find(map == unit01);
                    [row02, col02] = find(map == unit02);
                    dis = sqrt((row02-row01)^2 + (col02-col01)^2) * 0.25;
                    rtmp1 = [rtmp1, [rsc;dis]];
                end
            end
            
            if (rSig>0 && rSig<0.5)
                % Finding r_sc
                tunCur01 = tunMat(:,i);
                tunCur02 = tunMat(:,j);
                pref01 = find(tunCur01 == max(tunCur01));
                pref02 = find(tunCur02 == max(tunCur02));
                count01 = zeros(200,1);
                count02 = zeros(200,2);
                for trl = 1:200
                    tmp01 = S(pref01).trial(trl).counts;
                    tmp02 = S(pref02).trial(trl).counts;
                    count01(trl) = sum(tmp01(i,:));
                    count02(trl) = sum(tmp02(j,:));
                end
                % ====================================================!!!!!!!!!!!!!!!!!
                zsc = zscore([count01, count02]);
                rsc = corr(zsc);
                rsc = rsc(2);

                % finding distance
                unit01 = ch(i);
                unit02 = ch(j);
                if (unit01 ~= unit02)
                    [row01, col01] = find(map == unit01);
                    [row02, col02] = find(map == unit02);
                    dis = sqrt((row02-row01)^2 + (col02-col01)^2) * 0.25;
                    rtmp2 = [rtmp2, [rsc;dis]];
                end
            end
            
            if (rSig>-0.5 && rSig<0)
                % Finding r_sc
                tunCur01 = tunMat(:,i);
                tunCur02 = tunMat(:,j);
                pref01 = find(tunCur01 == max(tunCur01));
                pref02 = find(tunCur02 == max(tunCur02));
                count01 = zeros(200,1);
                count02 = zeros(200,2);
                for trl = 1:200
                    tmp01 = S(pref01).trial(trl).counts;
                    tmp02 = S(pref02).trial(trl).counts;
                    count01(trl) = sum(tmp01(i,:));
                    count02(trl) = sum(tmp02(j,:));
                end
                % ====================================================!!!!!!!!!!!!!!!!!
                zsc = zscore([count01, count02]);
                rsc = corr(zsc);
                rsc = rsc(2);
                
                % finding distance
                unit01 = ch(i);
                unit02 = ch(j);
                if (unit01 ~= unit02)
                    [row01, col01] = find(map == unit01);
                    [row02, col02] = find(map == unit02);
                    dis = sqrt((row02-row01)^2 + (col02-col01)^2) * 0.25;
                    rtmp3 = [rtmp3, [rsc;dis]];
                end
            end
            
            if (rSig<-0.5)
                % Finding r_sc
                tunCur01 = tunMat(:,i);
                tunCur02 = tunMat(:,j);
                pref01 = find(tunCur01 == max(tunCur01));
                pref02 = find(tunCur02 == max(tunCur02));
                count01 = zeros(200,1);
                count02 = zeros(200,2);
                for trl = 1:200
                    tmp01 = S(pref01).trial(trl).counts;
                    tmp02 = S(pref02).trial(trl).counts;
                    count01(trl) = sum(tmp01(i,:));
                    count02(trl) = sum(tmp02(j,:));
                end
                % ====================================================!!!!!!!!!!!!!!!!!
                zsc = zscore([count01, count02]);
                rsc = corr(zsc);
                rsc = rsc(2);

                % finding distance
                unit01 = ch(i);
                unit02 = ch(j);
                if (unit01 ~= unit02)
                    [row01, col01] = find(map == unit01);
                    [row02, col02] = find(map == unit02);
                    dis = sqrt((row02-row01)^2 + (col02-col01)^2) * 0.25;
                    rtmp4 = [rtmp4, [rsc;dis]];
                end
            end
        end
    end
    rtmp1(:,1) = [];
    rtmp2(:,1) = [];
    rtmp3(:,1) = [];
    rtmp4(:,1) = [];
    r1{imonkey} = rtmp1;
    r2{imonkey} = rtmp2;
    r3{imonkey} = rtmp3;
    r4{imonkey} = rtmp4;
    toc
end 

colors = {[.5 .5 .5],[1 .0 .0],[.0 .0 1],[.0 1 .0]};
thickness = [5,4,3,1.5];
leg = {'rSig > 0.5', '0 < rSig < 0.5', '-0.5 < rSig < 0', 'rSig < -0.5'};
for imonkey = 1:3
    rStruct = {r1{imonkey},r2{imonkey},r3{imonkey},r4{imonkey}};
    figure;
    for j = 1:4
        r = rStruct{j};
        x = unique(r(2,:));
        meanMat = [];
        stdMat = [];
        for i = 1:length(x)
            inds = find(r(2,:) == x(i));
            meanMat(i) = mean(r(1,inds));
            stdMat(i) = std(r(1,inds));
        end
        stdMat = stdMat/sqrt(200);
%         scatter(r(2,:),r(1,:));
        errorbar(x(1:3:end),meanMat(1:3:end),stdMat(1:3:end),'k','LineWidth',1,'LineStyle','none','HandleVisibility','off')
        hold on
        plot(x,meanMat,'color',colors{j},'LineWidth',thickness(j))
%         ylim([-0.1 0.3])
        xlabel("Distance Between Electrodes(mm)")
        ylabel("Spike Count Correlation(r_{sc})")
        title(sprintf("Monkey %d",imonkey))
    end
    legend(leg{1},leg{2},leg{3},leg{4});
end

figure;
for imonkey = 1:3
    r = r1{imonkey};
    x = unique(r(2,:));
    meanMat = [];
    stdMat = [];
    for i = 1:length(x)
        inds = find(r(2,:) == x(i));
        meanMat(i) = mean(r(1,inds));
        stdMat(i) = std(r(1,inds));
    end
    stdMat = stdMat/sqrt(200);
    subplot(1,3,imonkey)
    errorbar(x,meanMat,stdMat,'k','LineWidth',1,'LineStyle','none','HandleVisibility','off')
    plot(x,meanMat,'k','LineWidth',2)
    ylim([0 0.3])
    xlabel("Distance Between Electrodes(mm)")
    ylabel("Spike Count Correlation(r_{sc})")
    title(sprintf("Monkey %d (rSig > 0.5)",imonkey))
    grid on
end

%% Question 03 - Figure B&C =================================================
clear all; close all; clc;
binWidth = 20;
nMonkey = 3;
t = 0:binWidth:1000-binWidth;
gratsNum = 12;
monkeys = {'S_monkey1.mat','S_monkey2.mat','S_monkey3.mat'};
IND = {'IND_monkey1.mat','IND_monkey2.mat','IND_monkey3.mat'};
dataS = {'data_monkey1_gratings.mat','data_monkey2_gratings.mat','data_monkey3_gratings.mat'};

mat = {};
for imonkey = 1:nMonkey
    mat1 = [0;0;0];
    load(dataS{imonkey});
    load(IND{imonkey});
    load(monkeys{imonkey});

    map = data.MAP;
    removeInd = find(out == 0);
    ch = data.CHANNELS;
    ch(removeInd,:) = [];
    ch(:,2) = [];
    map(isnan(map)) = 0;

    % rSignal
    tunMat = zeros(gratsNum,length(find(out)));
    for i = 1:length(find(out))
        curve = zeros(gratsNum,1);
        for igrat = 1:gratsNum
            tmp = S(igrat).mean_FRs;
            curve(igrat) = mean(tmp(i,:));
        end
        tunMat(:,i) = curve;
    end
    tunCorr = corr(tunMat);

    % r_sc
    cntMat = zeros(200,length(find(out)));
    for i = 1:length(find(out))
        grat = find(tunMat(:,i) == max(tunMat(:,i)));
        cnt = zeros(200,1);
        for trl = 1:200
            sig = S(grat).trial(trl).counts;
            cnt(trl) = sum(sig(i,:));
        end
        cntMat(:,i) = cnt;
    end
    cntCorr = corr(cntMat);

    for i = 1:length(find(out))-1
        for j = i+1:length(find(out))
            unit01 = ch(i);
            unit02 = ch(j);
            if (unit01 ~= unit02)
                [row01, col01] = find(map == unit01);
                [row02, col02] = find(map == unit02);
                dis = 0.25 * sqrt((row02-row01)^2 + (col02-col01)^2);
            end
            rSig = tunCorr(i,j);
            rsc = cntCorr(i,j);
            mat1 = [mat1, [rSig;rsc;dis]]; % !!!!!!! This line is important !!!!!!!
        end
    end
    mat1(:,1) = [];
    mat{imonkey} = mat1;
end

% Making Figure B ============
for imonkey = 1:nMonkey
    matTmp = mat{imonkey};
    inds1 = find(matTmp(3,:) >= 0 & matTmp(3,:) <= 0.75);
    data1 = matTmp(:,inds1);
    x1 = unique(data1(1,:));
    meanMat1 = [];
    for i = 1:length(x1)
        inds = find(data1(1,:) == x1(i));
        meanMat1(i) = mean(data1(2,inds));
    end

    inds2 = find(matTmp(3,:) > 0.75 & matTmp(3,:) <= 1.25);
    data2 = matTmp(:,inds2);
    x2 = unique(data2(1,:));
    meanMat2 = [];
    for i = 1:length(x2)
        inds = find(data2(1,:) == x2(i));
        meanMat2(i) = mean(data2(2,inds));
    end

    inds3 = find(matTmp(3,:) > 1.25 & matTmp(3,:) <= 4);
    data3 = matTmp(:,inds3);
    x3 = unique(data3(1,:));
    meanMat3 = [];
    for i = 1:length(x3)
        inds = find(data3(1,:) == x3(i));
        meanMat3(i) = mean(data3(2,inds));
    end

    w = 100;
    meanMat1 = movmean(meanMat1,w);
    meanMat2 = movmean(meanMat2,w);
    meanMat3 = movmean(meanMat3,w);
    figure
    plot(x1,meanMat1,'k','LineWidth',3)
    hold on
    plot(x2,meanMat2,'r','LineWidth',2)
    hold on
    plot(x3,meanMat3,'b','LineWidth',1.5)
    xlabel("Orientation Tuning Similarity (rSig)")
    ylabel("Spike Count Correlation(r_{sc})")
    title(sprintf("Monkey %d",imonkey))
    grid on
    legend("d <= 0.75mm","0.75mm < d <= 1.25mm","1.25mm < d")
end

% Making Figure C ================
figure;
for imonkey = 1:nMonkey
    matTmp = mat{imonkey};
    dis = linspace(0.3,4,11);
    rSig = linspace(-0.9,0.9,11);
    rSig = flip(rSig);
    colorMat = zeros(length(rSig)-1,length(dis)-1);
    for i = 1:length(dis)-1
        dis01 = dis(i);
        dis02 = dis(i+1);
        for j = 1:length(rSig)-1
            rSig01 = rSig(j);
            rSig02 = rSig(j+1);
            inds = find(matTmp(3,:) >= dis01 & matTmp(3,:) <= dis02);
            tmp = matTmp(:,inds);
            inds = find(tmp(1,:) <= rSig01 & tmp(1,:) >= rSig02);
            tmp = tmp(:,inds);
            element = mean(tmp(2,:));
            colorMat(j,i) = element;
        end
    end
    subplot(1,3,imonkey)
    colormap(jet)
    imagesc(dis,rSig,colorMat);
    set(gca,'YDir','normal')
    colorbar
    xlabel("Distnace(mm)")
    ylabel("Orientation Tuning Similarity(rSig)")
    title(sprintf("Monkey %d",imonkey))
end


%% Question 5
clear all; close all; clc;






%% Functions
function RasterPlot(Trials, Fs)
    iTrials = size(Trials,1);
    Ts = 1/Fs * 1000; % in ms
    for i = 1:iTrials
        spike_times = find(Trials(i,:) == 1) *Ts;
        x = repmat(spike_times,3,1);
        y = nan(size(x));
        
        if ~isempty(y) % range from i-1 to i  
            y(1,:) = i - 1;
            y(2,:) = i;
        end
        
        plot(x, y, 'Color', 'k')
        hold on
    end
end






