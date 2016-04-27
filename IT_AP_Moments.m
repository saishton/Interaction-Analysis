function [IT_Struc,AP_Struc] = IT_AP_Moments(data,dir_ref)

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

%==Cut Extreme Data and Estimate Moments for Interaction Times==%
[F_times,X_times] = ecdf(times);
ccdf_times = 1-F_times;
Xrem = [X_times(end-2:end)];
X_times = X_times(2:end-3);
ccdf_times = ccdf_times(2:end-3);
dataMod = times(~ismember(times,Xrem));
times_test_data = sort(dataMod)';

M1 = mean(dataMod);
M2 = mean(dataMod.^2);
M3 = mean(dataMod.^3);

%==Estimate Parameters for Distributions on Interaction Times==%
times_lambda = M1;
times_ccdf_ex = expcdf(X_times,times_lambda,'upper');

%ML MoM does not work

[times_k,times_sigma,times_theta] = gpSolve(M1,M2,M3);
times_ccdf_gp = gpcdf(X_times,times_k,times_sigma,times_theta,'upper');

%WL MoM does not work

%==Extract GoF Data==%
grid = (min(times_test_data)-20):20:max(times_test_data);

times_z_ex = expcdf(grid,times_lambda);
times_z_gp = gpcdf(grid,times_k,times_sigma,times_theta);

times_stats_ex = testStatistics_d(times_test_data,times_z_ex,20);
times_stats_gp = testStatistics_d(times_test_data,times_z_gp,20);

%==Plotting==%
IntTimes_fig = figure();
hold on
plot(X_times,ccdf_times,'o')
plot(X_times,times_ccdf_ex)
plot(X_times,times_ccdf_gp)
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Contact Time (s)');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gen. Pareto','Location','southwest');
hold off
imagefilename = [dir_ref,'/InteractionTimes_Moments.png'];
print(imagefilename,'-dpng')
close(IntTimes_fig);

%==Cut Extreme Data and Estimate Moments for Activity Potential==%
[F_AP,X_AP] = ecdf(activityPot);
ap_ccdf = 1-F_AP;
Xrem = [];
X_AP = X_AP(2:end);
ap_ccdf = ap_ccdf(2:end);
dataMod = activityPot(~ismember(activityPot,Xrem));
ap_test_data = sort(dataMod)';

M1 = mean(dataMod);
M2 = mean(dataMod.^2);

%==Estimate Parameters for Distributions on Activity Potential==%
ap_lambda = M1;
ap_ccdf_ex = expcdf(X_AP,ap_lambda,'upper');

ap_b = (M2/M1)-M1;
ap_a = M1/ap_b;
ap_ccdf_gm = gamcdf(X_AP,ap_a,ap_b,'upper');

ap_sigma = M1*sqrt(2/pi);
ap_ccdf_rl = raylcdf(X_AP,ap_sigma,'upper');

ap_lnsig = sqrt(log(M2*(M1^-2)));
ap_lnmu = log(M1)-(0.5*ap_lnsig^2);
ap_ccdf_ln = logncdf(X_AP,ap_lnmu,ap_lnsig,'upper');

%==Extract GoF Data==%
ap_z_ex = expcdf(ap_test_data,ap_lambda);
ap_z_gm = gamcdf(ap_test_data,ap_a,ap_b);
ap_z_rl = raylcdf(ap_test_data,ap_sigma);
ap_z_ln = logncdf(ap_test_data,ap_lnmu,ap_lnsig);

ap_stats_ex = testStatistics(ap_test_data,ap_z_ex);
ap_stats_gm = testStatistics(ap_test_data,ap_z_gm);
ap_stats_rl = testStatistics(ap_test_data,ap_z_rl);
ap_stats_ln = testStatistics(ap_test_data,ap_z_ln);

ap_stats_ex.Root_MSE = sqrt(mean((ap_z_ex-ap_test_data).^2));
ap_stats_gm.Root_MSE = sqrt(mean((ap_z_gm-ap_test_data).^2));
ap_stats_rl.Root_MSE = sqrt(mean((ap_z_rl-ap_test_data).^2));
ap_stats_ln.Root_MSE = sqrt(mean((ap_z_ln-ap_test_data).^2));

TotSS = sum((ap_test_data-mean(ap_test_data)).^2);

ap_stats_ex.R_Squared = 1-sum((ap_test_data-ap_z_ex).^2)/TotSS;
ap_stats_gm.R_Squared = 1-sum((ap_test_data-ap_z_gm).^2)/TotSS;
ap_stats_rl.R_Squared = 1-sum((ap_test_data-ap_z_rl).^2)/TotSS;
ap_stats_ln.R_Squared = 1-sum((ap_test_data-ap_z_ln).^2)/TotSS;

%==Plotting==%
AP_fig = figure();
hold on
plot(X_AP,ap_ccdf,'o');
plot(X_AP,ap_ccdf_ex);
plot(X_AP,ap_ccdf_gm);
plot(X_AP,ap_ccdf_rl);
plot(X_AP,ap_ccdf_ln);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Activity Potential');
ylabel('CCDF');
axis([-inf,inf,1E-2,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
hold off
apfilename = [dir_ref,'/ActivityPotential_Moments.png'];
print(apfilename,'-dpng')
close(AP_fig);

%==Build and Return Relevant Data==%
times_struc_ex = struct('Scale',times_lambda);
times_struc_gp = struct('Shape',times_k,'Scale',times_sigma,'Location',times_theta);

ap_struc_ex = struct('Scale',ap_lambda);
ap_struc_gm = struct('Shape',ap_a,'Scale',ap_b);
ap_struc_rl = struct('Scale',ap_sigma);
ap_struc_ln = struct('Location',ap_lnmu,'Scale',ap_lnsig);

times_size = size(times,2);
ap_size = size(activityPot,2);

times_pvals_ex = pvals_ex_d(times_size,times_lambda,times_stats_ex,3,6,20);
times_pvals_gp = pvals_gp_d(times_size,times_k,times_sigma,times_theta,times_stats_gp,3,6,20);

ap_pvals_ex = pvals_ex(ap_size,ap_lambda,ap_stats_ex,0,6);
ap_pvals_gm = pvals_gm(ap_size,ap_a,ap_b,ap_stats_gm,0,6);
ap_pvals_rl = pvals_rl(ap_size,ap_sigma,ap_stats_rl,0,6);
ap_pvals_ln = pvals_ln(ap_size,ap_lnmu,ap_lnsig,ap_stats_ln,0,6);

times_EX = struct('Parameters',times_struc_ex,'Statistics',times_stats_ex,'pValues',times_pvals_ex);
times_GP = struct('Parameters',times_struc_gp,'Statistics',times_stats_gp,'pValues',times_pvals_gp);

AP_EX = struct('Parameters',ap_struc_ex,'Statistics',ap_stats_ex,'pValues',ap_pvals_ex);
AP_GM = struct('Parameters',ap_struc_gm,'Statistics',ap_stats_gm,'pValues',ap_pvals_gm);
AP_RL = struct('Parameters',ap_struc_rl,'Statistics',ap_stats_rl,'pValues',ap_pvals_rl);
AP_LN = struct('Parameters',ap_struc_ln,'Statistics',ap_stats_ln,'pValues',ap_pvals_ln);

IT_Struc = struct('Exponential',times_EX,'GenPareto',times_GP);
AP_Struc = struct('Exponential',AP_EX,'Gamma',AP_GM,'Rayleigh',AP_RL,'LogNormal',AP_LN);
end