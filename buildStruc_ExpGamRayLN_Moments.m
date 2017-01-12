function [Structure] = buildStruc_ExpGamRayLN_Moments(data,dir_ref,property_title,graph_title,cutExtreme,Ymin)

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
difference = diff(test_data);
difference = difference(difference>0);
res = min(difference);

%==Prepare MoM==%
M1 = mean(dataMod);
M2 = mean(dataMod.^2);

%==Extract parameters==%
ex_lambda = M1;
ccdf_ex = expcdf(X,ex_lambda,'upper');

gm_b = (M2/M1)-M1;
gm_a = M1/gm_b;
ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');

rl_sigma = M1*sqrt(2/pi);
ccdf_rl = raylcdf(X,rl_sigma,'upper');

ln_sigma = sqrt(log(M2*(M1^-2)));
ln_mu = log(M1)-(0.5*ln_sigma^2);
ccdf_ln = logncdf(X,ln_mu,ln_sigma,'upper');

%==Extract GoF data==%
z_ex = expcdf(test_data,ex_lambda);
z_gm = gamcdf(test_data,gm_a,gm_b);
z_rl = raylcdf(test_data,rl_sigma);
z_ln = logncdf(test_data,ln_mu,ln_sigma);

zp_ex = exppdf(test_data,ex_lambda);
zp_gm = gampdf(test_data,gm_a,gm_b);
zp_rl = raylpdf(test_data,rl_sigma);
zp_ln = lognpdf(test_data,ln_mu,ln_sigma);

stats_ex = testStatistics(test_data,z_ex,zp_ex,0);
stats_gm = testStatistics(test_data,z_gm,zp_gm,0);
stats_rl = testStatistics(test_data,z_rl,zp_rl,0);
stats_ln = testStatistics(test_data,z_ln,zp_ln,0);

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
imagefilename = [dir_ref,'/',property_title,'_Moments.png'];
print(imagefilename,'-dpng')
close(fig);

%==Build data structure==%
samplesize = max(size(data));

struc_ex = struct('Scale',ex_lambda);
struc_gm = struct('Shape',gm_a,'Scale',gm_b);
struc_rl = struct('Scale',rl_sigma);
struc_ln = struct('Location',ln_mu,'Scale',ln_sigma);

p_ex = pvals_ex(samplesize,ex_lambda,stats_ex,cutExtreme,MC_Power,res);
p_gm = pvals_gm(samplesize,gm_a,gm_b,stats_gm,cutExtreme,MC_Power,res);
p_rl = pvals_rl(samplesize,rl_sigma,stats_rl,cutExtreme,MC_Power,res);
p_ln = pvals_ln(samplesize,ln_mu,ln_sigma,stats_ln,cutExtreme,MC_Power,res);

EX = struct('Parameters',struc_ex,'Statistics',stats_ex,'pValues',p_ex);
GM = struct('Parameters',struc_gm,'Statistics',stats_gm,'pValues',p_gm);
RL = struct('Parameters',struc_rl,'Statistics',stats_rl,'pValues',p_rl);
LN = struct('Parameters',struc_ln,'Statistics',stats_ln,'pValues',p_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end