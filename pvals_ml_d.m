function [p_vals] = pvals_ml_d(dataLength,beta,gamma,Statistics,cut,n,step)

ModKol = Statistics.Modified_Kol;
ChiSq = Statistics.Chi_Squared;

num_MC = 10^n;

ModKol_stat = zeros(1,num_MC);
ChiSq_stat = zeros(1,num_MC);

parfor i=1:num_MC
	data = mlrnd(beta,gamma,dataLength,1);
    data = sort(data);
    if cut>0
        data(end-cut+1:end) = [];
    end 
    grid = min(data):step:max(data);
    CDF = ones(length(grid),1)-mlf(beta,1,-gamma*grid.^beta,6);
    thisfit = testStatistics_d(data,[0; CDF],step);
    
    ModKol_stat(i) = thisfit.Modified_Kol;
    ChiSq_stat(i) = thisfit.Chi_Squared;
end

[F_ModKol,X_ModKol] = ecdf(ModKol_stat);
[F_ChiSq,X_ChiSq] = ecdf(ChiSq_stat);

clearvars ModKol_stat ChiSq_stat

X_ModKol = [-1E99;X_ModKol(2:end);1E99];
X_ChiSq = [0;X_ChiSq(2:end);1E99];
F_ModKol = [0;F_ModKol(2:end);1];
F_ChiSq = [0;F_ChiSq(2:end);1];

p_ModKol = 1-interp1(X_ModKol,F_ModKol,ModKol,'next');
p_ChiSq = 1-interp1(X_ChiSq,F_ChiSq,ChiSq,'next');

p_vals = struct('Modified_Kol',p_ModKol,...
                'Chi_Squared',p_ChiSq);

end