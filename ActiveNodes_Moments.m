function [Structure] = ActiveNodes_Moments(data,dir_ref)

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

%==Cut Extreme Data and Estimate Moments==%
[F_Nodes,X_Nodes] = ecdf(numNodes);
ccdf_nodes = 1-F_Nodes;
Xrem = [X_Nodes(end-2:end)];
X_Nodes = X_Nodes(2:end-3);
ccdf_nodes = ccdf_nodes(2:end-3);
dataMod = numNodes(~ismember(numNodes,Xrem));
test_data = sort(dataMod)';

M1 = mean(dataMod);
M2 = mean(dataMod.^2);

%==Estimate Parameters for Distributions on Number of Active Nodes==%
nodes_lambda = M1;
nodes_ccdf_ex = expcdf(X_Nodes,nodes_lambda,'upper');

nodes_b = (M2/M1)-M1;
nodes_a = M1/nodes_b;
nodes_ccdf_gm = gamcdf(X_Nodes,nodes_a,nodes_b,'upper');

nodes_sigma = M1*sqrt(2/pi);
nodes_ccdf_rl = raylcdf(X_Nodes,nodes_sigma,'upper');

nodes_lnsig = sqrt(log(M2*(M1^-2)));
nodes_lnmu = log(M1)-(0.5*nodes_lnsig^2);
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

nodes_stats_ex.Root_MSE = sqrt(mean((nodes_z_ex-test_data).^2));
nodes_stats_gm.Root_MSE = sqrt(mean((nodes_z_gm-test_data).^2));
nodes_stats_rl.Root_MSE = sqrt(mean((nodes_z_rl-test_data).^2));
nodes_stats_ln.Root_MSE = sqrt(mean((nodes_z_ln-test_data).^2));

TotSS = sum((test_data-mean(test_data)).^2);

nodes_stats_ex.R_Squared = 1-sum((test_data-nodes_z_ex).^2)/TotSS;
nodes_stats_gm.R_Squared = 1-sum((test_data-nodes_z_gm).^2)/TotSS;
nodes_stats_rl.R_Squared = 1-sum((test_data-nodes_z_rl).^2)/TotSS;
nodes_stats_ln.R_Squared = 1-sum((test_data-nodes_z_ln).^2)/TotSS;

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
imagefilename = [dir_ref,'/ActiveNodes_Moments.png'];
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