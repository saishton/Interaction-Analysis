function [IT_Struc,AP_Struc] = IT_AP_FitTool(data,dir_ref)

number_rows = size(data,1);
number_people = max([data(:,2); data(:,3)]);
contact_time = 20;

%==Sort Data by ID==%
[~, order] = sort(data(:,3));
partsorteddata = data(order,:);
[~, order] = sort(partsorteddata(:,2));
sorteddata = partsorteddata(order,:);

%==Find Interaction Times==%
times = zeros(1,number_rows);
j = 1;
times_k = 1;
total_interactions = 0;
interactions = zeros(1,number_people);
step_vector = [contact_time 0 0];
while j<number_rows+1
    ID1 = sorteddata(j,2);
    ID2 = sorteddata(j,3);
    interactions(ID1) = interactions(ID1)+1;
    interactions(ID2) = interactions(ID2)+1;
    total_interactions = total_interactions+1;
    contact = contact_time;
    current_row = sorteddata(j,:);
    if j == number_rows
        next_row = [0 0 0]; 
    else
        next_row = sorteddata(j+1,:);
    end
    while isequal(next_row,current_row+step_vector)
        contact = contact+contact_time;
        j = j+1;
        current_row = sorteddata(j,:);
        if j == number_rows
            next_row = [0 0 0];
        else
            next_row = sorteddata(j+1,:);
        end
    end
    times(times_k) = contact;
    j = j+1;
    times_k = times_k+1;
end
times(times_k:end) = [];
interactions(interactions==0) = []; %Remove nodes where no interaction occurs
activityPot = interactions/total_interactions;

%==Fit Distributions for Interaction Times==%
[F_times,X_times] = ecdf(times);
ccdf_times = 1-F_times;
Xrem = [X_times(1);X_times(end-6:end)];
X_times = X_times(2:end-7);
ccdf_times = ccdf_times(2:end-7);

tau = mean(times); %Get rough estimate for exponential mean

fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[tau]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,gof_ex] = fit(X_times,ccdf_times,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[0.5 0.5]);
ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
[cf_ml,gof_ml] = fit(X_times,ccdf_times,ft_ml);
cv_ml = coeffvalues(cf_ml);

fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[0.5 0.5 0.5]);
ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
[cf_gp,gof_gp] = fit(X_times,ccdf_times,ft_gp);
cv_gp = coeffvalues(cf_gp);

fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[0.5 0.5]);
ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
[cf_wb,gof_wb] = fit(X_times,ccdf_times,ft_wb);
cv_wb = coeffvalues(cf_wb);

%==Extract Parameters==%
times_lambda = cv_ex(1);
times_ccdf_ex = expcdf(X_times,times_lambda,'upper');

times_beta = cv_ml(1);
times_gamma = cv_ml(2);
times_ccdf_ml = mlf(times_beta,1,-times_gamma*X_times.^times_beta,6);

times_k = cv_gp(1);
times_sigma = cv_gp(2);
times_theta = cv_gp(3);
times_ccdf_gp = gpcdf(X_times,times_k,times_sigma,times_theta,'upper');

times_a = cv_wb(1);
times_b = cv_wb(2);
times_ccdf_wb = wblcdf(X_times,times_a,times_b,'upper');

%==Extract GoF Data==%
dataMod = times(~ismember(times,Xrem));
test_data = sort(dataMod)';

times_z_ex = expcdf(test_data,times_lambda);
times_z_ml = ones(length(test_data),1)-mlf(times_beta,1,-times_gamma*test_data.^times_beta,6);
times_z_gp = gpcdf(test_data,times_k,times_sigma,times_theta);
times_z_wb = wblcdf(test_data,times_a,times_b);

times_stats_ex = testStatistics(test_data,times_z_ex);
times_stats_ml = testStatistics(test_data,times_z_ml);
times_stats_gp = testStatistics(test_data,times_z_gp);
times_stats_wb = testStatistics(test_data,times_z_wb);

times_stats_ex.Root_MSE = gof_ex.rmse;
times_stats_ml.Root_MSE = gof_ml.rmse;
times_stats_gp.Root_MSE = gof_gp.rmse;
times_stats_wb.Root_MSE = gof_wb.rmse;
times_stats_ex.R_Squared = gof_ex.rsquare;
times_stats_ml.R_Squared = gof_ml.rsquare;
times_stats_gp.R_Squared = gof_gp.rsquare;
times_stats_wb.R_Squared = gof_wb.rsquare;

%==Plotting==%
IntTimes_fig = figure();
hold on
plot(X_times,ccdf_times,'o')
plot(X_times,times_ccdf_ex)
plot(X_times,times_ccdf_ml)
plot(X_times,times_ccdf_gp)
plot(X_times,times_ccdf_wb)
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Contact Time (s)');
ylabel('CCDF');
axis([1E1,1E4,1E-5,1E0]);
legend('Data','Exponential','Mittag Leffler','Gen. Pareto','Weibull');
hold off
imagefilename = [dir_ref,'/InteractionTimes_FitTool.png'];
print(imagefilename,'-dpng')
close(IntTimes_fig);

%==Fit Distributions for Activity Potential==%
[F_AP,X_AP] = ecdf(activityPot);
ap_ccdf = 1-F_AP;
APXrem = [X_AP(1)];
X_AP = X_AP(2:end);
ap_ccdf = ap_ccdf(2:end);

fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
[cf_ex,gof_ex] = fit(X_AP,ap_ccdf,ft_ex);
cv_ex = coeffvalues(cf_ex);

fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
[cf_gm,gof_gm] = fit(X_AP,ap_ccdf,ft_gm);
cv_gm = coeffvalues(cf_gm);

fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
[cf_rl,gof_rl] = fit(X_AP,ap_ccdf,ft_rl);
cv_rl = coeffvalues(cf_rl);

fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[0 1]);
ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
[cf_ln,gof_ln] = fit(X_AP,ap_ccdf,ft_ln);
cv_ln = coeffvalues(cf_ln);

%==Extract Parameters==%
AP_lambda1 = cv_ex(1);
AP_ccdf_ex = expcdf(X_AP,AP_lambda1,'upper');

AP_a1 = cv_gm(1);
AP_b1 = cv_gm(2);
AP_ccdf_gm = gamcdf(X_AP,AP_a1,AP_b1,'upper');

AP_sigma1 = cv_rl(1);
AP_ccdf_rl = raylcdf(X_AP,AP_sigma1,'upper');

AP_lnmu = cv_ln(1);
AP_lnsig = cv_ln(2);
AP_ccdf_ln = logncdf(X_AP,AP_lnmu,AP_lnsig,'upper');

%==Extract GoF Data==%
dataMod = activityPot(~ismember(activityPot,APXrem));
test_data = sort(dataMod)';

AP_z_ex = expcdf(test_data,AP_lambda1);
AP_z_gm = gamcdf(test_data,AP_a1,AP_b1);
AP_z_rl = raylcdf(test_data,AP_sigma1);
AP_z_ln = logncdf(test_data,AP_lnmu,AP_lnsig);

AP_stats_ex = testStatistics(test_data,AP_z_ex);
AP_stats_gm = testStatistics(test_data,AP_z_gm);
AP_stats_rl = testStatistics(test_data,AP_z_rl);
AP_stats_ln = testStatistics(test_data,AP_z_ln);

AP_stats_ex.Root_MSE = gof_ex.rmse;
AP_stats_gm.Root_MSE = gof_gm.rmse;
AP_stats_rl.Root_MSE = gof_rl.rmse;
AP_stats_ln.Root_MSE = gof_ln.rmse;
AP_stats_ex.R_Squared = gof_ex.rsquare;
AP_stats_gm.R_Squared = gof_gm.rsquare;
AP_stats_rl.R_Squared = gof_rl.rsquare;
AP_stats_ln.R_Squared = gof_ln.rsquare;

%==Plotting==%
AP_fig = figure();
hold on
plot(X_AP,ap_ccdf,'o');
plot(X_AP,AP_ccdf_ex);
plot(X_AP,AP_ccdf_gm);
plot(X_AP,AP_ccdf_rl);
plot(X_AP,AP_ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Activity Potential');
ylabel('CCDF');
axis([1E-2,1E0,1E-2,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal');
hold off
apfilename = [dir_ref,'/ActivityPotential_FitTool.png'];
print(apfilename,'-dpng')
close(AP_fig);

%==Build and Return Relevant Data==%
times_struc_ex = struct('Scale',times_lambda);
times_struc_ml = struct('Stability',times_beta,'Scale',times_gamma);
times_struc_gp = struct('Shape',times_k,'Scale',times_sigma,'Location',times_theta);
times_struc_wb = struct('Scale',times_a,'Shape',times_b);

ap_struc_ex = struct('Scale',AP_lambda1);
ap_struc_gm = struct('Shape',AP_a1,'Scale',AP_b1);
ap_struc_rl = struct('Scale',AP_sigma1);
ap_struc_ln = struct('Location',AP_lnmu,'Scale',AP_lnsig);

times_EX = struct('Parameters',times_struc_ex,'Statistics',times_stats_ex);
times_ML = struct('Parameters',times_struc_ml,'Statistics',times_stats_ml);
times_GP = struct('Parameters',times_struc_gp,'Statistics',times_stats_gp);
times_WB = struct('Parameters',times_struc_wb,'Statistics',times_stats_wb);

AP_EX = struct('Parameters',ap_struc_ex,'Statistics',AP_stats_ex);
AP_GM = struct('Parameters',ap_struc_gm,'Statistics',AP_stats_gm);
AP_RL = struct('Parameters',ap_struc_rl,'Statistics',AP_stats_rl);
AP_LN = struct('Parameters',ap_struc_ln,'Statistics',AP_stats_ln);

IT_Struc = struct('Exponential',times_EX,'MittagLeffler',times_ML,'GenPareto',times_GP,'Weibull',times_WB);
AP_Struc = struct('Exponential',AP_EX,'Gamma',AP_GM,'Rayleigh',AP_RL,'LogNormal',AP_LN);
end