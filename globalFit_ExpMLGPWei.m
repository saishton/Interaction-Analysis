function [Structure] = globalFit_ExpMLGPWei(allData,mins,maxs,numSplits)

cutExtreme = 3;
MC_Power = 6;
MC_Power_ML = 3;

fieldNames = fieldnames(allData);
numData = length(fieldNames);

for p=1:numData
    thisname = fieldNames{p};
    thisdata = allData.(thisname);
    dataSize.(thisname) = length(thisdata);
    [~,X] = ecdf(thisdata);
    if cutExtreme>0
        Xrem = [X(end+1-cutExtreme:end)];
    else
        Xrem = [];
    end
    dataMod = thisdata(~ismember(thisdata,Xrem));
    test_data.(thisname) = sort(dataMod)';
end

EX1 = linspace(mins(1),maxs(1),numSplits);
ML1 = linspace(mins(2),maxs(2),numSplits);
ML2 = linspace(mins(3),maxs(3),numSplits);
GP1 = linspace(mins(4),maxs(4),numSplits);
GP2 = linspace(mins(5),maxs(5),numSplits);
GP3 = linspace(mins(6),maxs(6),numSplits);
WB1 = linspace(mins(7),maxs(7),numSplits);
WB2 = linspace(mins(8),maxs(8),numSplits);

EX_BestStats = Inf(1,5);
EX_Stats = zeros(5,5);
EX1_BestCoords = zeros(1,5);
for p=1:numSplits
    thisEX1 = EX1(p);
    EX1_thisCoord = thisEX1*ones(1,5);
    thisStats = zeros(numData,5);
    parfor j=1:numData
        thisName = fieldNames{j};
        thisData = test_data.(thisName);
        z = expcdf(thisData,thisEX1);
        thisTest = testStatistics(thisData,z);
        thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
    end
    if numData > 1
        thisSum = sum(thisStats);
    else
        thisSum = thisStats;
    end
    EX1_BestCoords(thisSum<EX_BestStats) = EX1_thisCoord(thisSum<EX_BestStats);
    fullSum = ones(5,1)*thisSum;
    EX_Stats(thisSum<EX_BestStats,:) = fullSum(thisSum<EX_BestStats,:);
    EX_BestStats = min(EX_BestStats,thisSum);
end

ML_BestStats = Inf(1,5);
ML_Stats = zeros(5,5);
ML1_BestCoords = zeros(1,5);
ML2_BestCoords = zeros(1,5);
for p=1:numSplits
    thisML1 = ML1(p);
    for q=1:numSplits
        thisML2 = ML2(q);
        ML1_thisCoord = thisML1*ones(1,5);
        ML2_thisCoord = thisML2*ones(1,5);
        thisStats = zeros(numData,5);
        parfor j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = ones(length(thisData),1)-mlf(thisML1,1,-thisML2*thisData.^thisML1,6);
            thisTest = testStatistics(thisData,z);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        ML1_BestCoords(thisSum<ML_BestStats) = ML1_thisCoord(thisSum<ML_BestStats);
        ML2_BestCoords(thisSum<ML_BestStats) = ML2_thisCoord(thisSum<ML_BestStats);
        fullSum = ones(5,1)*thisSum;
        ML_Stats(thisSum<ML_BestStats,:) = fullSum(thisSum<ML_BestStats,:);
        ML_BestStats = min(ML_BestStats,thisSum);
    end
end

GP_BestStats = Inf(1,5);
GP_Stats = zeros(5,5);
GP1_BestCoords = zeros(1,5);
GP2_BestCoords = zeros(1,5);
GP3_BestCoords = zeros(1,5);
for p=1:numSplits
    thisGP1 = GP1(p);
    for q=1:numSplits
        thisGP2 = GP2(q);
        for r = 1:numSplits
            thisGP3 = GP3(r);
            GP1_thisCoord = thisGP1*ones(1,5);
            GP2_thisCoord = thisGP2*ones(1,5);
            GP3_thisCoord = thisGP3*ones(1,5);
            thisStats = zeros(numData,5);
            parfor j=1:numData
                thisName = fieldNames{j};
                thisData = test_data.(thisName);
                z = raylcdf(thisData,thisGP1);
                thisTest = testStatistics(thisData,z);
                thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
            end
            if numData > 1
                thisSum = sum(thisStats);
            else
                thisSum = thisStats;
            end
            GP1_BestCoords(thisSum<GP_BestStats) = GP1_thisCoord(thisSum<GP_BestStats);
            GP2_BestCoords(thisSum<GP_BestStats) = GP2_thisCoord(thisSum<GP_BestStats);
            GP3_BestCoords(thisSum<GP_BestStats) = GP3_thisCoord(thisSum<GP_BestStats);
            fullSum = ones(5,1)*thisSum;
            GP_Stats(thisSum<GP_BestStats,:) = fullSum(thisSum<GP_BestStats,:);
            GP_BestStats = min(GP_BestStats,thisSum);
        end
    end
end

WB_BestStats = Inf(1,5);
WB_Stats = zeros(5,5);
WB1_BestCoords = zeros(1,5);
WB2_BestCoords = zeros(1,5);
for p=1:numSplits
    thisWB1 = WB1(p);
    for q=1:numSplits
        thisWB2 = WB2(q);
        WB1_thisCoord = thisWB1*ones(1,5);
        WB2_thisCoord = thisWB2*ones(1,5);
        thisStats = zeros(numData,5);
        parfor j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = wblcdf(thisData,thisWB1,thisWB2);
            thisTest = testStatistics(thisData,z);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        WB1_BestCoords(thisSum<WB_BestStats) = WB1_thisCoord(thisSum<WB_BestStats);
        WB2_BestCoords(thisSum<WB_BestStats) = WB2_thisCoord(thisSum<WB_BestStats);
        fullSum = ones(5,1)*thisSum;
        WB_Stats(thisSum<WB_BestStats,:) = fullSum(thisSum<WB_BestStats,:);
        WB_BestStats = min(WB_BestStats,thisSum);
    end
end

statsMatrix = [EX_BestStats;ML_BestStats;GP_BestStats;WB_BestStats];
[~,idx] = min(statsMatrix);
structureNames = {'KolD_Pri';'CvM_Pri';'Kuiper_Pri';'Watson_Pri';'AD_Pri'};

for i=1:5
    thisName = structureNames{i};
    if idx(i)==1
        Structure.(thisName).Distribution = struct('Type','Exponential',...
                                            'Scale',EX1_BestCoords(i)...
                                            );
        stat = EX_Stats(i,:);
    elseif idx(i)==2
        Structure.(thisName).Distribution = struct('Type','Mittag-Leffler',...
                                            'Stability',ML1_BestCoords(i),...
                                            'Scale',ML2_BestCoords(i)...
                                            );
        stat = ML_Stats(i,:);
    elseif idx(i)==3
        Structure.(thisName).Distribution = struct('Type','Generalised Pareto',...
                                            'Shape',GP1_BestCoords(i),...
                                            'Scale',GP2_BestCoords(i),...
                                            'Location',GP3_BestCoords(i)...
                                            );
        stat = GP_Stats(i,:);
    elseif idx(i)==4
        Structure.(thisName).Distribution = struct('Type','Weibull',...
                                            'Scale',WB1_BestCoords(i),...
                                            'Shape',WB2_BestCoords(i)...
                                            );
        stat = WB_Stats(i,:);
    end
    Structure.(thisName).Statistics.Kolmogorov_D = stat(1);
    Structure.(thisName).Statistics.Cramer_von_Mises = stat(2);
    Structure.(thisName).Statistics.Kuiper = stat(3);
    Structure.(thisName).Statistics.Watson = stat(4);
    Structure.(thisName).Statistics.Anderson_Darling = stat(5);
end

for  i=1:5
    thisName = structureNames{i};
    thisStructure = Structure.(thisName);
    thisDistribution = thisStructure.Distribution;
    thisStatistics = thisStructure.Statistics;
    if strcmp(thisDistribution.Type,'Exponential')
        Structure.(thisName).pValues = global_pvals_ex(dataSize,thisDistribution.Scale,thisStatistics,cutExtreme,MC_Power);
    elseif strcmp(thisDistribution.Type,'Mittag-Leffler')
        Structure.(thisName).pValues = global_pvals_ml(dataSize,thisDistribution.Stability,thisDistribution.Scale,thisStatistics,cutExtreme,MC_Power_ML);
    elseif strcmp(thisDistribution.Type,'Generalised Pareto')
        Structure.(thisName).pValues = global_pvals_gp(dataSize,thisDistribution.Shape,thisDistribution.Scale,thisDistribution.Location,thisStatistics,cutExtreme,MC_Power);
    elseif strcmp(thisDistribution.Type,'Weibull')
        Structure.(thisName).pValues = global_pvals_wb(dataSize,thisDistribution.Scale,thisDistribution.Shape,thisStatistics,cutExtreme,MC_Power);
    end
end