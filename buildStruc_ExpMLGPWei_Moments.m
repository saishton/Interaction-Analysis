function [Structure] = buildStruc_ExpMLGPWei_Moments(data,dir_ref,property_title,graph_title,cutExtreme,Ymin)

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

%==Prepare MoM==%
M1 = mean(dataMod);
M2 = mean(dataMod.^2);
M3 = mean(dataMod.^3);

%==Extract parameters==%
ex_lambda = M1;
ccdf_ex = expcdf(X,ex_lambda,'upper');

[gp_k,gp_sigma,gp_theta] = gpSolve(M1,M2,M3);
ccdf_gp = gpcdf(X,gp_k,gp_sigma,gp_theta,'upper');

%==Extract GoF data==%
z_ex = expcdf(test_data,ex_lambda);
z_gp = gpcdf(test_data,gp_k,gp_sigma,gp_theta);

zp_ex = exppdf(test_data,ex_lambda);
zp_gp = gppdf(test_data,gp_k,gp_sigma,gp_theta);

stats_ex = testStatistics(test_data,z_ex,zp_ex);
stats_gp = testStatistics(test_data,z_gp,zp_gp);

%==Plotting==%
fig = figure();
hold on
plot(X,ccdf,'o');
plot(X,ccdf_ex);
plot(X,ccdf_gp);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(graph_title);
ylabel('CCDF');
axis([-inf,inf,Ymin,1E0]);
legend('Data','Exponential','Gen. Pareto','Location','southwest');
imagefilename = [dir_ref,'/',property_title,'_Moments.png'];
print(imagefilename,'-dpng')
close(fig);

%==Build data structure==%
samplesize = max(size(data));

struc_ex = struct('Scale',ex_lambda);
struc_gp = struct('Shape',gp_k,'Scale',gp_sigma,'Location',gp_theta);

p_ex = pvals_ex(samplesize,ex_lambda,stats_ex,cutExtreme,MC_Power);
p_gp = pvals_gp(samplesize,gp_k,gp_sigma,gp_theta,stats_gp,cutExtreme,MC_Power);

EX = struct('Parameters',struc_ex,'Statistics',stats_ex,'pValues',p_ex);
GP = struct('Parameters',struc_gp,'Statistics',stats_gp,'pValues',p_gp);

Structure = struct('Exponential',EX,'GenPareto',GP);
end