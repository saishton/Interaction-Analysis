function [Structure] = buildStruc_ExpGamRayLN_FitTool(data,dir_ref,property_title,graph_title,cutExtreme)

MC_Power = 6;

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

ex_lambda_start = M1;
gm_b_start = (M2/M1)-M1;
gm_a_start = M1/gm_b_start;
rl_sigma_start = M1*sqrt(2/pi);
ln_sigma_start = sqrt(log(M2*(M1^-2)));
ln_mu_start = log(M1)-(0.5*ln_sigma_start^2);

%==Fit distributions==%
fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,~] = fit(X,ccdf,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[gm_a_start gm_b_start]);
ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
[cf_gm,~] = fit(X,ccdf,ft_gm);
cv_gm = coeffvalues(cf_gm);

fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[rl_sigma_start]);
ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
[cf_rl,~] = fit(X,ccdf,ft_rl);
cv_rl = coeffvalues(cf_rl);

fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[ln_mu_start ln_sigma_start]);
ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
[cf_ln,~] = fit(X,ccdf,ft_ln);
cv_ln = coeffvalues(cf_ln);

%==Extract parameters==%
ex_lambda = cv_ex(1);
ccdf_ex = expcdf(X,ex_lambda,'upper');

gm_a = cv_gm(1);
gm_b = cv_gm(2);
ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');

rl_sigma = cv_rl(1);
ccdf_rl = raylcdf(X,rl_sigma,'upper');

ln_mu = cv_ln(1);
ln_sigma = cv_ln(2);
ccdf_ln = logncdf(X,ln_mu,ln_sigma,'upper');

%==Extract GoF data==%
z_ex = expcdf(test_data,ex_lambda);
z_gm = gamcdf(test_data,gm_a,gm_b);
z_rl = raylcdf(test_data,rl_sigma);
z_ln = logncdf(test_data,ln_mu,ln_sigma);

stats_ex = testStatistics(test_data,z_ex);
stats_gm = testStatistics(test_data,z_gm);
stats_rl = testStatistics(test_data,z_rl);
stats_ln = testStatistics(test_data,z_ln);

%==Plotting==%
fig = figure();
hold on
plot(X,ccdf,'o');
plot(X,ccdf_ex);
plot(X,ccdf_gm);
plot(X,ccdf_rl);
plot(X,ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(graph_title);
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
imagefilename = [dir_ref,'/',property_title,'_FitTool.png'];
print(imagefilename,'-dpng')
close(fig);

%==Build data structure==%
samplesize = max(size(data));

struc_ex = struct('Scale',ex_lambda);
struc_gm = struct('Shape',gm_a,'Scale',gm_b);
struc_rl = struct('Scale',rl_sigma);
struc_ln = struct('Location',ln_mu,'Scale',ln_sigma);

p_ex = pvals_ex(samplesize,ex_lambda,stats_ex,cutExtreme,MC_Power);
p_gm = pvals_gm(samplesize,gm_a,gm_b,stats_gm,cutExtreme,MC_Power);
p_rl = pvals_rl(samplesize,rl_sigma,stats_rl,cutExtreme,MC_Power);
p_ln = pvals_ln(samplesize,ln_mu,ln_sigma,stats_ln,cutExtreme,MC_Power);

EX = struct('Parameters',struc_ex,'Statistics',stats_ex,'pValues',p_ex);
GM = struct('Parameters',struc_gm,'Statistics',stats_gm,'pValues',p_gm);
RL = struct('Parameters',struc_rl,'Statistics',stats_rl,'pValues',p_rl);
LN = struct('Parameters',struc_ln,'Statistics',stats_ln,'pValues',p_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end