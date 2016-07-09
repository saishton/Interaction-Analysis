function [Structure] = NumComp_MLE(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

CompCount = zeros(1,num_times);

parfor m=1:num_times
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
[thisCompCount,~,~] = networkComponents(thisadj);
CompCount(m) = thisCompCount;
end

%==Fit Distributions for Number of Components==%
[F_comp,X_comp] = ecdf(CompCount);
ccdf_comp = 1-F_comp;
Xrem = [X_comp(end-2:end)];
X_comp = X_comp(2:end-3);
ccdf_comp = ccdf_comp(2:end-3);
dataMod = CompCount(~ismember(CompCount,Xrem));
test_data = sort(dataMod)';

mod_test_data = test_data;
mod_test_data(mod_test_data==0)=1E-99;

phat_ex = mle(mod_test_data,'distribution','Exponential');
phat_gm = mle(mod_test_data,'distribution','Gamma');
phat_rl = mle(mod_test_data,'distribution','Rayleigh');
phat_ln = mle(mod_test_data,'distribution','Lognormal');

%==Extract Parameters==%
comp_lambda = phat_ex(1);
comp_ccdf_ex = expcdf(X_comp,comp_lambda,'upper');

comp_a = phat_gm(1);
comp_b = phat_gm(2);
comp_ccdf_gm = gamcdf(X_comp,comp_a,comp_b,'upper');

comp_sigma = phat_rl(1);
comp_ccdf_rl = raylcdf(X_comp,comp_sigma,'upper');

comp_lnmu = phat_ln(1);
comp_lnsig = phat_ln(2);
comp_ccdf_ln = logncdf(X_comp,comp_lnmu,comp_lnsig,'upper');

%==Extract GoF Data==%
comp_z_ex = expcdf(test_data,comp_lambda);
comp_z_gm = gamcdf(test_data,comp_a,comp_b);
comp_z_rl = raylcdf(test_data,comp_sigma);
comp_z_ln = logncdf(test_data,comp_lnmu,comp_lnsig);

comp_stats_ex = testStatistics(test_data,comp_z_ex);
comp_stats_gm = testStatistics(test_data,comp_z_gm);
comp_stats_rl = testStatistics(test_data,comp_z_rl);
comp_stats_ln = testStatistics(test_data,comp_z_ln);

%==Plotting==%
numCompfig = figure();
hold on
plot(X_comp,ccdf_comp,'o');
plot(X_comp,comp_ccdf_ex);
plot(X_comp,comp_ccdf_gm);
plot(X_comp,comp_ccdf_rl);
plot(X_comp,comp_ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of Components');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
imagefilename = [dir_ref,'/NumberComponents_MLE.png'];
print(imagefilename,'-dpng')
close(numCompfig);

%==Build and Return Relevant Data==%
comp_struc_ex = struct('Scale',comp_lambda);
comp_struc_gm = struct('Shape',comp_a,'Scale',comp_b);
comp_struc_rl = struct('Scale',comp_sigma);
comp_struc_ln = struct('Location',comp_lnmu,'Scale',comp_lnsig);

comp_pvals_ex = pvals_ex(num_times,comp_lambda,comp_stats_ex,3,6);
comp_pvals_gm = pvals_gm(num_times,comp_a,comp_b,comp_stats_gm,3,6);
comp_pvals_rl = pvals_rl(num_times,comp_sigma,comp_stats_rl,3,6);
comp_pvals_ln = pvals_ln(num_times,comp_lnmu,comp_lnsig,comp_stats_ln,3,6);

EX = struct('Parameters',comp_struc_ex,'Statistics',comp_stats_ex,'pValues',comp_pvals_ex);
GM = struct('Parameters',comp_struc_gm,'Statistics',comp_stats_gm,'pValues',comp_pvals_gm);
RL = struct('Parameters',comp_struc_rl,'Statistics',comp_stats_rl,'pValues',comp_pvals_rl);
LN = struct('Parameters',comp_struc_ln,'Statistics',comp_stats_ln,'pValues',comp_pvals_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end