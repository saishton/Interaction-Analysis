function [p_vals] = pvals_wb(dataLength,a,b,Statistics,cut,n)

KolDPlus = Statistics.Kolmogorov_D_Plus;
KolDMinus = Statistics.Kolmogorov_D_Minus;
KolD = Statistics.Kolmogorov_D;
CvM = Statistics.Cramer_von_Mises;
Kuiper = Statistics.Kuiper;
Watson = Statistics.Watson;
AD = Statistics.Anderson_Darling;

num_MC = 10^n;

KolDPlus_stat = zeros(1,num_MC);
KolDMinus_stat = zeros(1,num_MC);
KolD_stat = zeros(1,num_MC);
CvM_stat = zeros(1,num_MC);
Kuiper_stat = zeros(1,num_MC);
Watson_stat = zeros(1,num_MC);
AD_stat = zeros(1,num_MC);

parfor i=1:num_MC
    data = wblrnd(a,b,dataLength,1);
    data = sort(data);
    if cut>0
        data(end-cut+1:end) = [];
    end 
    CDF = wblcdf(data,a,b);
    thisfit = testStatistics(data,CDF);
    
    KolDPlus_stat(i) = thisfit.Kolmogorov_D_Plus;
    KolDMinus_stat(i) = thisfit.Kolmogorov_D_Minus;
    KolD_stat(i) = thisfit.Kolmogorov_D;
    CvM_stat(i) = thisfit.Cramer_von_Mises;
    Kuiper_stat(i) = thisfit.Kuiper;
    Watson_stat(i) = thisfit.Watson;
    AD_stat(i) = thisfit.Anderson_Darling;
end

[F_KolDPlus,X_KolDPlus] = ecdf(KolDPlus_stat);
[F_KolDMinus,X_KolDMinus] = ecdf(KolDMinus_stat);
[F_KolD,X_KolD] = ecdf(KolD_stat);
[F_CvM,X_CvM] = ecdf(CvM_stat);
[F_Kuiper,X_Kuiper] = ecdf(Kuiper_stat);
[F_Watson,X_Watson] = ecdf(Watson_stat);
[F_AD,X_AD] = ecdf(AD_stat);

clearvars KolDPlus_stat KolDMinus_stat KolD_stat CvM_stat Kuiper_stat Watson_stat AD_stat

p_KolDPlus = 1-interp1(X_KolDPlus(2:end),F_KolDPlus(2:end),KolDPlus,'next');
p_KolDMinus = 1-interp1(X_KolDMinus(2:end),F_KolDMinus(2:end),KolDMinus,'next');
p_KolD = 1-interp1(X_KolD(2:end),F_KolD(2:end),KolD,'next');
p_CvM = 1-interp1(X_CvM(2:end),F_CvM(2:end),CvM,'next');
p_Kuiper = 1-interp1(X_Kuiper(2:end),F_Kuiper(2:end),Kuiper,'next');
p_Watson = 1-interp1(X_Watson(2:end),F_Watson(2:end),Watson,'next');
p_AD = 1-interp1(X_AD(2:end),F_AD(2:end),AD,'next');

p_vals = struct('Kolmogorov_D_Plus',p_KolDPlus,...
                'Kolmogorov_D_Minus',p_KolDMinus,...
                'Kolmogorov_D',p_KolD,...
                'Cramer_von_Mises',p_CvM,...
                'Kuiper',p_Kuiper,...
                'Watson',p_Watson,...
                'Anderson_Darling',p_AD);

end