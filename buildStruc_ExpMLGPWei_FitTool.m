function [Structure] = buildStruc_ExpMLGPWei_FitTool(data,dir_ref,property_title,graph_title,cutExtreme)

MC_Power = 6;
MC_Power_ML = 3;

%==Prepare data==%
[F,X] = ecdf(data);
ccdf = 1-F;
if cutExtreme>0
    Xrem = [X(end+1-cutExtreme:end)];
else
    Xrem = [];
end
X = X(2:end-cutExtreme);
ccdf = ccdf(2:end-cutExtreme);
dataMod = data(~ismember(data,Xrem));
test_data = sort(dataMod)';

%==Choose initial values (using MoM)==%
M1 = mean(dataMod);
M2 = mean(dataMod.^2);
M3 = mean(dataMod.^3);

ex_lambda_start = M1;
ml_beta_start = 0.5;
ml_gamma_start = 0.5;
[gp_k_start,gp_sigma_start,gp_theta_start] = gpSolve(M1,M2,M3);
wb_a_start = 0.5;
wb_b_start = 0.5;

%==Fit distributions==%
fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,~] = fit(X,ccdf,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[ml_beta_start ml_gamma_start]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
[cf_ml,~] = fit(X,ccdf,ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[gp_k_start gp_sigma_start gp_theta_start]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
[cf_gp,~] = fit(X,ccdf,ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[wb_a_start wb_b_start]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
[cf_wb,~] = fit(X,ccdf,ft_wb);
cv_wb = coeffvalues(cf_wb);

%==Extract parameters==%
ex_lambda = cv_ex(1);
ccdf_ex = expcdf(X,ex_lambda,'upper');

ml_beta = cv_ml(1);
ml_gamma = cv_ml(2);
ccdf_ml = mlf(ml_beta,1,-ml_gamma*X.^ml_beta,6);

gp_k = cv_gp(1);
gp_sigma = cv_gp(2);
gp_theta = cv_gp(3);
ccdf_gp = gpcdf(X,gp_k,gp_sigma,gp_theta,'upper');

wb_a = cv_wb(1);
wb_b = cv_wb(2);
ccdf_wb = wblcdf(X,wb_a,wb_b,'upper');

%==Extract GoF data==%
z_ex = expcdf(test_data,ex_lambda);
z_ml = ones(length(test_data),1)-mlf(ml_beta,1,-ml_gamma*test_data.^ml_beta,6);
z_gp = gpcdf(test_data,gp_k,gp_sigma,gp_theta);
z_wb = wblcdf(test_data,wb_a,wb_b);

stats_ex = testStatistics(test_data,z_ex);
stats_ml = testStatistics(test_data,z_ml);
stats_gp = testStatistics(test_data,z_gp);
stats_wb = testStatistics(test_data,z_wb);

%==Plotting==%
fig = figure();
hold on
plot(X,ccdf,'o');
plot(X,ccdf_ex);
plot(X,ccdf_ml);
plot(X,ccdf_gp);
plot(X,ccdf_wb);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(graph_title);
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Mittag Leffler','Gen. Pareto','Weibull','Location','southwest');
imagefilename = [dir_ref,'/',property_title,'_FitTool.png'];
print(imagefilename,'-dpng')
close(fig);

%==Build data structure==%
samplesize = max(size(data));

struc_ex = struct('Scale',ex_lambda);
struc_ml = struct('Stability',ml_beta,'Scale',ml_gamma);
struc_gp = struct('Shape',gp_k,'Scale',gp_sigma,'Location',gp_theta);
struc_wb = struct('Scale',wb_a,'Shape',wb_b);

p_ex = pvals_ex(samplesize,ex_lambda,stats_ex,cutExtreme,MC_Power);
p_ml = pvals_ml(samplesize,ml_beta,ml_gamma,stats_ml,cutExtreme,MC_Power_ML);
p_gp = pvals_gp(samplesize,gp_k,gp_sigma,gp_theta,stats_gp,cutExtreme,MC_Power);
p_wb = pvals_wb(samplesize,wb_a,wb_b,stats_wb,cutExtreme,MC_Power);

EX = struct('Parameters',struc_ex,'Statistics',stats_ex,'pValues',p_ex);
ML = struct('Parameters',struc_ml,'Statistics',stats_ml,'pValues',p_ml);
GP = struct('Parameters',struc_gp,'Statistics',stats_gp,'pValues',p_gp);
WB = struct('Parameters',struc_wb,'Statistics',stats_wb,'pValues',p_wb);

Structure = struct('Exponential',EX,'MittagLeffler',ML,'GenPareto',GP,'Weibull',WB);
end