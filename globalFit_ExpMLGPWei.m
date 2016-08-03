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

ML_AllStats = Inf(numSplits*numSplits,7);
for p=1:numSplits
    thisML1 = ML1(p);
    parfor q=1:numSplits
        thisML2 = ML2(q);
        thisStats = zeros(numData,7);
        for j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = ones(length(thisData),1)-mlf(thisML1,1,-thisML2*thisData.^thisML1,6);
            zprime = 
            thisTest = testStatistics(thisData,z,zprime);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        pML_AllStats(q,:) = thisSum;
    end
    low = (p-1)*numSplits+1;
    high = p*numSplits;
    ML_AllStats(low:high,:) = pML_AllStats;
end
[ML_BestStats,idx] = min(ML_AllStats);
ML_Stats = zeros(7,7);
ML1_BestCoords = zeros(1,7);
ML2_BestCoords = zeros(1,7);
parfor i=1:7
    ML_Stats(i,:) = ML_AllStats(idx(i),:);
    p = floor(idx(i)/numSplits)+1;
    q = rem(idx(i),numSplits);
    if q==0
        p = p-1;
        q = numSplits;
    end
    ML1_BestCoords(i) = ML1(p);
    ML2_BestCoords(i) = ML2(q);
end

GP_AllStats = Inf(numSplits^3,7);
for p=1:numSplits
    thisGP1 = GP1(p);
    for q=1:numSplits
        thisGP2 = GP2(q);
        parfor r = 1:numSplits
            thisGP3 = GP3(r);
            thisStats = zeros(numData,7);
            for j=1:numData
                thisName = fieldNames{j};
                thisData = test_data.(thisName);
                z = gpcdf(thisData,thisGP1);
                zprime = gppdf(thisData,thisGP1);
                thisTest = testStatistics(thisData,z,zprime);
                thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
            end
            if numData > 1
                thisSum = sum(thisStats);
            else
                thisSum = thisStats;
            end
            rGP_AllStats(r,:) = thisSum;
        end
        low = (q-1)*numSplits+1;
        high = q*numSplits;
        qGP_AllStats(low:high,:) = rGP_AllStats;
    end
    low = (p-1)*(numSplits^2)+1;
    high = p*(numSplits^2);
    GP_AllStats(low:high,:) = qGP_AllStats;
end
[GP_BestStats,idx] = min(GP_AllStats);
GP_Stats = zeros(7,7);
GP1_BestCoords = zeros(1,7);
GP2_BestCoords = zeros(1,7);
GP3_BestCoords = zeros(1,7);
parfor i=1:7
    GP_Stats(i,:) = GP_AllStats(idx(i),:);
    p = floor(idx(i)/(numSplits^2))+1;
    qr = rem(idx(i),(numSplits^2));
    if qr==0
        p = p-1;
        qr = numSplits^2;
    end
    q = floor(qr/numSplits)+1;
    r = rem(qr,numSplits);
    if r==0
        q = p-1;
        r = numSplits;
    end
    GP1_BestCoords(i) = GP1(p);
    GP2_BestCoords(i) = GP2(q);
    GP3_BestCoords(i) = GP3(r);
end

WB_AllStats = Inf(numSplits*numSplits,7);
for p=1:numSplits
    thisWB1 = WB1(p);
    parfor q=1:numSplits
        thisWB2 = WB2(q);
        thisStats = zeros(numData,7);
        for j=1:numData
            thisName = fieldNames{j};
            thisData = test_data.(thisName);
            z = wblcdf(thisData,thisWB1,thisWB2);
            zprime = wblpdf(thisData,thisWB1,thisWB2);
            thisTest = testStatistics(thisData,z,zprime);
            thisStats(j,:) = [thisTest.Kolmogorov_D,thisTest.Cramer_von_Mises,thisTest.Kuiper,thisTest.Watson,thisTest.Anderson_Darling,thisTest.Kullback_Leibler,thisTest.Jensen_Shannon];
        end
        if numData > 1
            thisSum = sum(thisStats);
        else
            thisSum = thisStats;
        end
        pWB_AllStats(q,:) = thisSum;
    end
    low = (p-1)*numSplits+1;
    high = p*numSplits;
    WB_AllStats(low:high,:) = pWB_AllStats;
end
[WB_BestStats,idx] = min(WB_AllStats);
WB_Stats = zeros(7,7);
WB1_BestCoords = zeros(1,7);
WB2_BestCoords = zeros(1,7);
parfor i=1:7
    WB_Stats(i,:) = WB_AllStats(idx(i),:);
    p = floor(idx(i)/numSplits)+1;
    q = rem(idx(i),numSplits);
    if q==0
        p = p-1;
        q = numSplits;
    end
    WB1_BestCoords(i) = WB1(p);
    WB2_BestCoords(i) = WB2(q);
end

statsMatrix = [EX_BestStats;ML_BestStats;GP_BestStats;WB_BestStats];
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
    elseif strcmp(thisDistribution.Type,'Mittag-Leffler')
        Structure.(thisName).pValues = global_pvals_ml(dataSize,thisDistribution.Stability,thisDistribution.Scale,thisStatistics,cutExtreme,MC_Power_ML);
    elseif strcmp(thisDistribution.Type,'Generalised Pareto')
        Structure.(thisName).pValues = global_pvals_gp(dataSize,thisDistribution.Shape,thisDistribution.Scale,thisDistribution.Location,thisStatistics,cutExtreme,MC_Power);
    elseif strcmp(thisDistribution.Type,'Weibull')
        Structure.(thisName).pValues = global_pvals_wb(dataSize,thisDistribution.Scale,thisDistribution.Shape,thisStatistics,cutExtreme,MC_Power);
    end
end