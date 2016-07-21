function [p_vals] = global_pvals_ml(dataLength,beta,gamma,Statistics,cut,n)

KolD = Statistics.Kolmogorov_D;
CvM = Statistics.Cramer_von_Mises;
Kuiper = Statistics.Kuiper;
Watson = Statistics.Watson;
AD = min(Statistics.Anderson_Darling,1E99);

num_MC = 10^n;
fieldNames = fieldnames(dataLength);
numData = length(fieldNames);

KolD_stat = zeros(1,num_MC);
CvM_stat = zeros(1,num_MC);
Kuiper_stat = zeros(1,num_MC);
Watson_stat = zeros(1,num_MC);
AD_stat = zeros(1,num_MC);

parfor i=1:num_MC
    thisStats = zeros(numData,5);
    for j=1:numData
        thisName = fieldNames{j};
        thisLength = dataLength.(thisName);
        if thisLength >0
            data = mlrnd(beta,gamma,thisLength,1);
            data = sort(data);
            if cut>0
                data(end-cut+1:end) = [];
            end    
            CDF = ones(length(data),1)-mlf(beta,1,-gamma*data.^beta,6);
            thisfit = testStatistics(data,CDF);
            thisStats(j,:) = [thisfit.Kolmogorov_D,thisfit.Cramer_von_Mises,thisfit.Kuiper,thisfit.Watson,thisfit.Anderson_Darling];
        end
    end
    if numData > 1
        thisSum = sum(thisStats);
    else
        thisSum = thisStats;
    end
    KolD_stat(i) = thisSum(1);
    CvM_stat(i) = thisSum(2);
    Kuiper_stat(i) = thisSum(3);
    Watson_stat(i) = thisSum(4);
    AD_stat(i) = thisSum(5);
end

[F_KolD,X_KolD] = ecdf(KolD_stat);
[F_CvM,X_CvM] = ecdf(CvM_stat);
[F_Kuiper,X_Kuiper] = ecdf(Kuiper_stat);
[F_Watson,X_Watson] = ecdf(Watson_stat);
[F_AD,X_AD] = ecdf(AD_stat);

clearvars KolD_stat CvM_stat Kuiper_stat Watson_stat AD_stat

X_KolD = [-1E99;X_KolD(2:end);1E99];
X_CvM = [-1E99;X_CvM(2:end);1E99];
X_Kuiper = [-1E99;X_Kuiper(2:end);1E99];
X_Watson = [-1E99;X_Watson(2:end);1E99];
X_AD = [-1E99;X_AD(2:end);1E99];

F_KolD = [0;F_KolD(2:end);1];
F_CvM = [0;F_CvM(2:end);1];
F_Kuiper = [0;F_Kuiper(2:end);1];
F_Watson = [0;F_Watson(2:end);1];
F_AD = [0;F_AD(2:end);1];

p_KolD = 1-interp1(X_KolD,F_KolD,KolD,'next');
p_CvM = 1-interp1(X_CvM,F_CvM,CvM,'next');
p_Kuiper = 1-interp1(X_Kuiper,F_Kuiper,Kuiper,'next');
p_Watson = 1-interp1(X_Watson,F_Watson,Watson,'next');
p_AD = 1-interp1(X_AD,F_AD,AD,'next');

p_vals = struct('Kolmogorov_D',p_KolD,...
                'Cramer_von_Mises',p_CvM,...
                'Kuiper',p_Kuiper,...
                'Watson',p_Watson,...
                'Anderson_Darling',p_AD);

end