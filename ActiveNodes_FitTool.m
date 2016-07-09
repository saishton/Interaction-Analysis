function [Structure] = ActiveNodes_FitTool(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

numNodes = zeros(1,num_times);

for m=1:num_times
    thisactive = zeros(1,num_people);
    current_time = (m-1)*contact_time;
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            thisactive(person1) = 1;
            thisactive(person2) = 1;
        end
    end
    numNodes(m) = sum(thisactive);
end

%==Fit Distributions for Number of Active Nodes==%
[F_Nodes,X_Nodes] = ecdf(numNodes);
ccdf_nodes = 1-F_Nodes;
Xrem = [X_Nodes(end-2:end)];
X_Nodes = X_Nodes(2:end-3);
ccdf_nodes = ccdf_nodes(2:end-3);
dataMod = numNodes(~ismember(numNodes,Xrem));
test_data = sort(dataMod)';

M1 = mean(dataMod);
M2 = mean(dataMod.^2);
b_start = (M2/M1)-M1;
a_start = M1/b_start;
lnsig_start = sqrt(log(M2*(M1^-2)));
lnmu_start = log(M1)-(0.5*lnsig_start^2);

nodes_fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
nodes_ft_ex = fittype('expcdf(x,lambda,''upper'')','options',nodes_fo_ex);
[nodes_cf_ex,nodes_gof_ex] = fit(X_Nodes,ccdf_nodes,nodes_ft_ex);
nodes_cv_ex = coeffvalues(nodes_cf_ex);

nodes_fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[a_start b_start]);
nodes_ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',nodes_fo_gm);
[nodes_cf_gm,nodes_gof_gm] = fit(X_Nodes,ccdf_nodes,nodes_ft_gm);
nodes_cv_gm = coeffvalues(nodes_cf_gm);

nodes_fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
nodes_ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',nodes_fo_rl);
[nodes_cf_rl,nodes_gof_rl] = fit(X_Nodes,ccdf_nodes,nodes_ft_rl);
nodes_cv_rl = coeffvalues(nodes_cf_rl);

nodes_fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[lnmu_start lnsig_start]);
nodes_ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',nodes_fo_ln);
[nodes_cf_ln,nodes_gof_ln] = fit(X_Nodes,ccdf_nodes,nodes_ft_ln);
nodes_cv_ln = coeffvalues(nodes_cf_ln);

%==Extract Parameters==%
nodes_lambda = nodes_cv_ex(1);
nodes_ccdf_ex = expcdf(X_Nodes,nodes_lambda,'upper');

nodes_a = nodes_cv_gm(1);
nodes_b = nodes_cv_gm(2);
nodes_ccdf_gm = gamcdf(X_Nodes,nodes_a,nodes_b,'upper');

nodes_sigma = nodes_cv_rl(1);
nodes_ccdf_rl = raylcdf(X_Nodes,nodes_sigma,'upper');

nodes_lnmu = nodes_cv_ln(1);
nodes_lnsig = nodes_cv_ln(2);
nodes_ccdf_ln = logncdf(X_Nodes,nodes_lnmu,nodes_lnsig,'upper');

%==Extract GoF Data==%
nodes_z_ex = expcdf(test_data,nodes_lambda);
nodes_z_gm = gamcdf(test_data,nodes_a,nodes_b);
nodes_z_rl = raylcdf(test_data,nodes_sigma);
nodes_z_ln = logncdf(test_data,nodes_lnmu,nodes_lnsig);

nodes_stats_ex = testStatistics(test_data,nodes_z_ex);
nodes_stats_gm = testStatistics(test_data,nodes_z_gm);
nodes_stats_rl = testStatistics(test_data,nodes_z_rl);
nodes_stats_ln = testStatistics(test_data,nodes_z_ln);

nodes_stats_ex.Root_MSE = nodes_gof_ex.rmse;
nodes_stats_gm.Root_MSE = nodes_gof_gm.rmse;
nodes_stats_rl.Root_MSE = nodes_gof_rl.rmse;
nodes_stats_ln.Root_MSE = nodes_gof_ln.rmse;
nodes_stats_ex.R_Squared = nodes_gof_ex.rsquare;
nodes_stats_gm.R_Squared = nodes_gof_gm.rsquare;
nodes_stats_rl.R_Squared = nodes_gof_rl.rsquare;
nodes_stats_ln.R_Squared = nodes_gof_ln.rsquare;

%==Plotting==%
nodesactivefig = figure();
hold on
plot(X_Nodes,ccdf_nodes,'o');
plot(X_Nodes,nodes_ccdf_ex);
plot(X_Nodes,nodes_ccdf_gm);
plot(X_Nodes,nodes_ccdf_rl);
plot(X_Nodes,nodes_ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of Nodes Active');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
imagefilename = [dir_ref,'/ActiveNodes_FitTool.png'];
print(imagefilename,'-dpng')
close(nodesactivefig);

%==Build and Return Relevant Data==%
nodes_struc_ex = struct('Scale',nodes_lambda);
nodes_struc_gm = struct('Shape',nodes_a,'Scale',nodes_b);
nodes_struc_rl = struct('Scale',nodes_sigma);
nodes_struc_ln = struct('Location',nodes_lnmu,'Scale',nodes_lnsig);

nodes_pvals_ex = pvals_ex(num_times,nodes_lambda,nodes_stats_ex,3,6);
nodes_pvals_gm = pvals_gm(num_times,nodes_a,nodes_b,nodes_stats_gm,3,6);
nodes_pvals_rl = pvals_rl(num_times,nodes_sigma,nodes_stats_rl,3,6);
nodes_pvals_ln = pvals_ln(num_times,nodes_lnmu,nodes_lnsig,nodes_stats_ln,3,6);

EX = struct('Parameters',nodes_struc_ex,'Statistics',nodes_stats_ex,'pValues',nodes_pvals_ex);
GM = struct('Parameters',nodes_struc_gm,'Statistics',nodes_stats_gm,'pValues',nodes_pvals_gm);
RL = struct('Parameters',nodes_struc_rl,'Statistics',nodes_stats_rl,'pValues',nodes_pvals_rl);
LN = struct('Parameters',nodes_struc_ln,'Statistics',nodes_stats_ln,'pValues',nodes_pvals_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end