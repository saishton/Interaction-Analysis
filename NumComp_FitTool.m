function [Structure] = NumComp_FitTool(data,dir_ref)

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

M1 = mean(dataMod);
M2 = mean(dataMod.^2);
b_start = (M2/M1)-M1;
a_start = M1/b_start;
lnsig_start = sqrt(log(M2*(M1^-2)));
lnmu_start = log(M1)-(0.5*lnsig_start^2);

comp_fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
comp_ft_ex = fittype('expcdf(x,lambda,''upper'')','options',comp_fo_ex);
[comp_cf_ex,comp_gof_ex] = fit(X_comp,ccdf_comp,comp_ft_ex);
comp_cv_ex = coeffvalues(comp_cf_ex);

comp_fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[a_start b_start]);
comp_ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',comp_fo_gm);
[comp_cf_gm,comp_gof_gm] = fit(X_comp,ccdf_comp,comp_ft_gm);
comp_cv_gm = coeffvalues(comp_cf_gm);

comp_fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
comp_ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',comp_fo_rl);
[comp_cf_rl,comp_gof_rl] = fit(X_comp,ccdf_comp,comp_ft_rl);
comp_cv_rl = coeffvalues(comp_cf_rl);

comp_fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[lnmu_start lnsig_start]);
comp_ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',comp_fo_ln);
[comp_cf_ln,comp_gof_ln] = fit(X_comp,ccdf_comp,comp_ft_ln);
comp_cv_ln = coeffvalues(comp_cf_ln);

%==Extract Parameters==%
comp_lambda = comp_cv_ex(1);
comp_ccdf_ex = expcdf(X_comp,comp_lambda,'upper');

comp_a = comp_cv_gm(1);
comp_b = comp_cv_gm(2);
comp_ccdf_gm = gamcdf(X_comp,comp_a,comp_b,'upper');

comp_sigma = comp_cv_rl(1);
comp_ccdf_rl = raylcdf(X_comp,comp_sigma,'upper');

comp_lnmu = comp_cv_ln(1);
comp_lnsig = comp_cv_ln(2);
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

comp_stats_ex.Root_MSE = comp_gof_ex.rmse;
comp_stats_gm.Root_MSE = comp_gof_gm.rmse;
comp_stats_rl.Root_MSE = comp_gof_rl.rmse;
comp_stats_ln.Root_MSE = comp_gof_ln.rmse;
comp_stats_ex.R_Squared = comp_gof_ex.rsquare;
comp_stats_gm.R_Squared = comp_gof_gm.rsquare;
comp_stats_rl.R_Squared = comp_gof_rl.rsquare;
comp_stats_ln.R_Squared = comp_gof_ln.rsquare;

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
imagefilename = [dir_ref,'/NumberComponents_FitTool.png'];
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