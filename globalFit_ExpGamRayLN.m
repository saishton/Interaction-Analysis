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

EX_AllStats = Inf(numSplits,7);
parfor p=1:numSplits
    thisEX1 = EX1(p);
    thisStats = zeros(numData,7);
    for j=1:numData
        thisName = fieldNames{j};
        thisData = test_data.(thisName);
        z = expcdf(thisData,thisEX1);
        zprime = exppdf(thisData,thisEX1);
        thisTest = testStatistics(thisData,z,zprime);
        thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
    end
    if numData > 1
        thisSum = sum(thisStats);
    else
        thisSum = thisStats;
    end
    EX_AllStats(p,:) = thisSum;
end
[EX_BestStats,idx] = min(EX_AllStats);
EX_Stats = zeros(7,7);
EX1_BestCoords = zeros(1,7);
parfor i=1:7
    EX_Stats(i,:) = EX_AllStats(idx(i),:);
    EX1_BestCoords(i) = EX1(idx(i));
end

GM_AllStats = Inf(numSplits*numSplits,7);
for p=1:numSplits
    thisGM1 = GM1(p);
    parfor q=1:numSplits
        thisGM2 = GM2(q);
        thisStats = zeros(numData,7);
        for j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = gamcdf(thisData,thisGM1,thisGM2);
            zprime = gampdf(thisData,thisGM1,thisGM2);
            thisTest = testStatistics(thisData,z,zprime);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        pGM_AllStats(q,:) = thisSum;
    end
    low = (p-1)*numSplits+1;
    high = p*numSplits;
    GM_AllStats(low:high,:) = pGM_AllStats;
end
[GM_BestStats,idx] = min(GM_AllStats);
GM_Stats = zeros(7,7);
GM1_BestCoords = zeros(1,7);
GM2_BestCoords = zeros(1,7);
parfor i=1:7
    GM_Stats(i,:) = GM_AllStats(idx(i),:);
    p = floor(idx(i)/numSplits)+1;
    q = rem(idx(i),numSplits);
    if q==0
        p = p-1;
        q = numSplits;
    end
    GM1_BestCoords(i) = GM1(p);
    GM2_BestCoords(i) = GM2(q);
end

RL_AllStats = Inf(numSplits,7);
parfor p=1:numSplits
    thisRL1 = RL1(p);
    thisStats = zeros(numData,7);
    for j=1:numData
        thisName = fieldNames{j};
        thisData = test_data.(thisName);
        z = raylcdf(thisData,thisRL1);
        zprime = raylpdf(thisData,thisRL1);
        thisTest = testStatistics(thisData,z,zprime);
        thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
    end
    if numData > 1
        thisSum = sum(thisStats);
    else
        thisSum = thisStats;
    end
    RL_AllStats(p,:) = thisSum;
end
[RL_BestStats,idx] = min(RL_AllStats);
RL_Stats = zeros(7,7);
RL1_BestCoords = zeros(1,7);
parfor i=1:7
    RL_Stats(i,:) = RL_AllStats(idx(i),:);
    RL1_BestCoords(i) = RL1(idx(i));
end

LN_AllStats = Inf(numSplits*numSplits,7);
for p=1:numSplits
    thisLN1 = LN1(p);
    parfor q=1:numSplits
        thisLN2 = LN2(q);
        thisStats = zeros(numData,7);
        for j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = logncdf(thisData,thisLN1,thisLN2);
            zprime = lognpdf(thisData,thisLN1,thisLN2);
            thisTest = testStatistics(thisData,z,zprime);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        pLN_AllStats(q,:) = thisSum;
    end
    low = (p-1)*numSplits+1;
    high = p*numSplits;
    LN_AllStats(low:high,:) = pLN_AllStats;
end
[LN_BestStats,idx] = min(LN_AllStats);
LN_Stats = zeros(7,7);
LN1_BestCoords = zeros(1,7);
LN2_BestCoords = zeros(1,7);
parfor i=1:7
    LN_Stats(i,:) = LN_AllStats(idx(i),:);
    p = floor(idx(i)/numSplits)+1;
    q = rem(idx(i),numSplits);
    if q==0
        p = p-1;
        q = numSplits;
    end
    LN1_BestCoords(i) = LN1(p);
    LN2_BestCoords(i) = LN2(q);
end

statsMatrix = [EX_BestStats;GM_BestStats;RL_BestStats;LN_BestStats];
[~,idx] = min(statsMatrix);
structureNames = {'KolD_Pri';'CvM_Pri';'Kuiper_Pri';'Watson_Pri';'AD_Pri';'KL_Pri';'JS_Pri';};

for i=1:7
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
    Structure.(thisName).Statistics.Kullback_Leibler = stat(6);
    Structure.(thisName).Statistics.Jensen_Shannon = stat(7);
end

for  i=1:7
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