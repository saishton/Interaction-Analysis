function [Structure] = buildStruc_ExpGamRayLN_MLE(data,dir_ref,property_title,graph_title,cutExtreme,Ymin)

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

%==Perform MLEs==%
mod_test_data = test_data;
mod_test_data(mod_test_data==0)=[];

phat_ex = mle(mod_test_data,'distribution','Exponential');
phat_gm = mle(mod_test_data,'distribution','Gamma');
phat_rl = mle(mod_test_data,'distribution','Rayleigh');
phat_ln = mle(mod_test_data,'distribution','Lognormal');

%==Extract parameters==%
ex_lambda = phat_ex(1);
ccdf_ex = expcdf(X,ex_lambda,'upper');

gm_a = phat_gm(1);
gm_b = phat_gm(2);
ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');

rl_sigma = phat_rl(1);
ccdf_rl = raylcdf(X,rl_sigma,'upper');

ln_mu = phat_ln(1);
ln_sigma = phat_ln(2);
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
axis([-inf,inf,Ymin,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
imagefilename = [dir_ref,'/',property_title,'_MLE.png'];
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