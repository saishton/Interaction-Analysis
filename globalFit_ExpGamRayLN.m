function [Structure] = globalFit_ExpGamRayLN(allData,mins,maxs,numSplits)

cutExtreme = 3;
MC_Power = 6;

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
GM1 = linspace(mins(2),maxs(2),numSplits);
GM2 = linspace(mins(3),maxs(3),numSplits);
RL1 = linspace(mins(4),maxs(4),numSplits);
LN1 = linspace(mins(5),maxs(5),numSplits);
LN2 = linspace(mins(6),maxs(6),numSplits);

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

GM_BestStats = Inf(1,5);
GM_Stats = zeros(5,5);
GM1_BestCoords = zeros(1,5);
GM2_BestCoords = zeros(1,5);
for p=1:numSplits
    thisGM1 = GM1(p);
    for q=1:numSplits
        thisGM2 = GM2(q);
        GM1_thisCoord = thisGM1*ones(1,5);
        GM2_thisCoord = thisGM2*ones(1,5);
        thisStats = zeros(numData,5);
        parfor j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = gamcdf(thisData,thisGM1,thisGM2);
            thisTest = testStatistics(thisData,z);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        GM1_BestCoords(thisSum<GM_BestStats) = GM1_thisCoord(thisSum<GM_BestStats);
        GM2_BestCoords(thisSum<GM_BestStats) = GM2_thisCoord(thisSum<GM_BestStats);
        fullSum = ones(5,1)*thisSum;
        GM_Stats(thisSum<GM_BestStats,:) = fullSum(thisSum<GM_BestStats,:);
        GM_BestStats = min(GM_BestStats,thisSum);
    end
end

RL_BestStats = Inf(1,5);
RL_Stats = zeros(5,5);
RL1_BestCoords = zeros(1,5);
for p=1:numSplits
    thisRL1 = RL1(p);
    RL1_thisCoord = thisRL1*ones(1,5);
    thisStats = zeros(numData,5);
    parfor j=1:numData
        thisName = fieldNames{j};
        thisData = test_data.(thisName);
        z = raylcdf(thisData,thisRL1);
        thisTest = testStatistics(thisData,z);
        thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
    end
    if numData > 1
        thisSum = sum(thisStats);
    else
        thisSum = thisStats;
    end
    RL1_BestCoords(thisSum<RL_BestStats) = RL1_thisCoord(thisSum<RL_BestStats);
    fullSum = ones(5,1)*thisSum;
    RL_Stats(thisSum<RL_BestStats,:) = fullSum(thisSum<RL_BestStats,:);
    RL_BestStats = min(RL_BestStats,thisSum);
end

LN_BestStats = Inf(1,5);
LN_Stats = zeros(5,5);
LN1_BestCoords = zeros(1,5);
LN2_BestCoords = zeros(1,5);
for p=1:numSplits
    thisLN1 = LN1(p);
    for q=1:numSplits
        thisLN2 = LN2(q);
        LN1_thisCoord = thisLN1*ones(1,5);
        LN2_thisCoord = thisLN2*ones(1,5);
        thisStats = zeros(numData,5);
        parfor j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = gamcdf(thisData,thisLN1,thisLN2);
            thisTest = testStatistics(thisData,z);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        LN1_BestCoords(thisSum<LN_BestStats) = LN1_thisCoord(thisSum<LN_BestStats);
        LN2_BestCoords(thisSum<LN_BestStats) = LN2_thisCoord(thisSum<LN_BestStats);
        fullSum = ones(5,1)*thisSum;
        LN_Stats(thisSum<LN_BestStats,:) = fullSum(thisSum<LN_BestStats,:);
        LN_BestStats = min(LN_BestStats,thisSum);
    end
end

statsMatrix = [EX_BestStats;GM_BestStats;RL_BestStats;LN_BestStats];
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
        Structure.(thisName).Distribution = struct('Type','Gamma',...
                                            'Shape',GM1_BestCoords(i),...
                                            'Scale',GM2_BestCoords(i)...
                                            );
        stat = GM_Stats(i,:);
    elseif idx(i)==3
        Structure.(thisName).Distribution = struct('Type','Rayleigh',...
                                            'Shape',RL1_BestCoords(i)...
                                            );
        stat = RL_Stats(i,:);
    elseif idx(i)==4
        Structure.(thisName).Distribution = struct('Type','Log-Normal',...
                                            'Shape',LN1_BestCoords(i),...
                                            'Scale',LN2_BestCoords(i)...
                                            );
        stat = LN_Stats(i,:);
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
    elseif strcmp(thisDistribution.Type,'Gamma')
        Structure.(thisName).pValues = global_pvals_gm(dataSize,thisDistribution.Shape,thisDistribution.Scale,thisStatistics,cutExtreme,MC_Power);
    elseif strcmp(thisDistribution.Type,'Rayleigh')
        Structure.(thisName).pValues = global_pvals_rl(dataSize,thisDistribution.Shape,thisStatistics,cutExtreme,MC_Power);
    elseif strcmp(thisDistribution.Type,'Log-Normal')
        Structure.(thisName).pValues = global_pvals_ln(dataSize,thisDistribution.Shape,thisDistribution.Scale,thisStatistics,cutExtreme,MC_Power);
    end
end