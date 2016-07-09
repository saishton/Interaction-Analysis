function [Structure] = AL_MLE(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

clustering = zeros(1,num_times);
numlinks = zeros(1,num_times);

for m=1:num_times
    thisadj = zeros(num_people);
    current_time = (m-1)*contact_time;
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            thisadj(person1,person2) = 1;
            thisadj(person2,person1) = 1;
        end
    end
    adj2 = thisadj^2;
    adj3 = thisadj^3;
    adj2sum = sum(sum(adj2));
    contrip = adj2sum - trace(adj2);
    if contrip==0
        clustering(m) = 0;
    else
        clustering(m) = trace(adj3)/contrip;
    end
    adjsum = sum(sum(thisadj));
    numlinks(m) = adjsum/2;
end

%==Fit Distributions for Number of Active Links==%
[F_links,X_links] = ecdf(numlinks);
ccdf_links = 1-F_links;
Xrem = [X_links(end-2:end)];
X_links = X_links(2:end-3);
ccdf_links = ccdf_links(2:end-3);
dataMod = numlinks(~ismember(numlinks,Xrem));
links_test_data = sort(dataMod)';

mod_test_data = links_test_data;
mod_test_data(mod_test_data==0)=1E-99;

phat_ex = mle(mod_test_data,'distribution','Exponential');
phat_gm = mle(mod_test_data,'distribution','Gamma');
phat_rl = mle(mod_test_data,'distribution','Rayleigh');
phat_ln = mle(mod_test_data,'distribution','Lognormal');

%==Extract Parameters==%
links_lambda = phat_ex(1);
links_ccdf_ex = expcdf(X_links,links_lambda,'upper');

links_a = phat_gm(1);
links_b = phat_gm(2);
links_ccdf_gm = gamcdf(X_links,links_a,links_b,'upper');

links_sigma = phat_rl(1);
links_ccdf_rl = raylcdf(X_links,links_sigma,'upper');

links_lnmu = phat_ln(1);
links_lnsig = phat_ln(2);
links_ccdf_ln = logncdf(X_links,links_lnmu,links_lnsig,'upper');

%==Extract GoF Data==%
links_z_ex = expcdf(links_test_data,links_lambda);
links_z_gm = gamcdf(links_test_data,links_a,links_b);
links_z_rl = raylcdf(links_test_data,links_sigma);
links_z_ln = logncdf(links_test_data,links_lnmu,links_lnsig);

links_stats_ex = testStatistics(links_test_data,links_z_ex);
links_stats_gm = testStatistics(links_test_data,links_z_gm);
links_stats_rl = testStatistics(links_test_data,links_z_rl);
links_stats_ln = testStatistics(links_test_data,links_z_ln);

%==Plotting==%
linksactivefig = figure();
hold on
plot(X_links,ccdf_links,'o');
plot(X_links,links_ccdf_ex);
plot(X_links,links_ccdf_gm);
plot(X_links,links_ccdf_rl);
plot(X_links,links_ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of Links Active');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
imagefilename = [dir_ref,'/LinkActivations_MLE.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Build and Return Relevant Data==%
links_struc_ex = struct('Scale',links_lambda);
links_struc_gm = struct('Shape',links_a,'Scale',links_b);
links_struc_rl = struct('Scale',links_sigma);
links_struc_ln = struct('Location',links_lnmu,'Scale',links_lnsig);

links_pvals_ex = pvals_ex(num_times,links_lambda,links_stats_ex,3,6);
links_pvals_gm = pvals_gm(num_times,links_a,links_b,links_stats_gm,3,6);
links_pvals_rl = pvals_rl(num_times,links_sigma,links_stats_rl,3,6);
links_pvals_ln = pvals_ln(num_times,links_lnmu,links_lnsig,links_stats_ln,3,6);

EX = struct('Parameters',links_struc_ex,'Statistics',links_stats_ex,'pValues',links_pvals_ex);
GM = struct('Parameters',links_struc_gm,'Statistics',links_stats_gm,'pValues',links_pvals_gm);
RL = struct('Parameters',links_struc_rl,'Statistics',links_stats_rl,'pValues',links_pvals_rl);
LN = struct('Parameters',links_struc_ln,'Statistics',links_stats_ln,'pValues',links_pvals_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end