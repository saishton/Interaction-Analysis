function [IT_Struc,AP_Struc] = IT_AP_MLE(data,dir_ref)

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
Xrem = [X_times(end-2:end)];
X_times = X_times(2:end-3);
ccdf_times = ccdf_times(2:end-3);
dataMod = times(~ismember(times,Xrem));
test_data = sort(dataMod)';

ml_pdf = @(mod_test_data,bet,gam)(-mlf(bet,0,-gam*mod_test_data.^bet,6)./mod_test_data);
ml_start = [1E0 1E-1];

mod_test_data = test_data;
mod_test_data(mod_test_data==0)=1E-99;

phat_ex = mle(mod_test_data,'distribution','Exponential');
phat_wb = mle(mod_test_data,'distribution','Weibull');

%==Extract Parameters==%
times_lambda = phat_ex(1);
times_ccdf_ex = expcdf(X_times,times_lambda,'upper');


times_a = phat_wb(1);
times_b = phat_wb(2);
times_ccdf_wb = wblcdf(X_times,times_a,times_b,'upper');

%==Extract GoF Data==%
times_z_ex = expcdf(test_data,times_lambda);
times_z_wb = wblcdf(test_data,times_a,times_b);

times_stats_ex = testStatistics(test_data,times_z_ex);
times_stats_wb = testStatistics(test_data,times_z_wb);

%==Plotting==%
IntTimes_fig = figure();
hold on
plot(X_times,ccdf_times,'o')
plot(X_times,times_ccdf_ex)
plot(X_times,times_ccdf_wb)
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Contact Time (s)');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Weibull','Location','southwest');
hold off
imagefilename = [dir_ref,'/InteractionTimes_MLE.png'];
print(imagefilename,'-dpng')
close(IntTimes_fig);

%==Fit Distributions for Activity Potential==%
[F_AP,X_AP] = ecdf(activityPot);
ap_ccdf = 1-F_AP;
APXrem = [X_AP(1)];
X_AP = X_AP(2:end);
ap_ccdf = ap_ccdf(2:end);
dataMod2 = activityPot(~ismember(activityPot,APXrem));
test_data2 = sort(dataMod2)';

mod_test_data2 = test_data2;
mod_test_data2(mod_test_data2==0)=1E-99;

phat_ex2 = mle(mod_test_data2,'distribution','Exponential');
phat_gm2 = mle(mod_test_data2,'distribution','Gamma');
phat_rl2 = mle(mod_test_data2,'distribution','Rayleigh');
phat_ln2 = mle(mod_test_data2,'distribution','Lognormal');

%==Extract Parameters==%
AP_lambda1 = phat_ex2(1);
AP_ccdf_ex = expcdf(X_AP,AP_lambda1,'upper');

AP_a1 = phat_gm2(1);
AP_b1 = phat_gm2(2);
AP_ccdf_gm = gamcdf(X_AP,AP_a1,AP_b1,'upper');

AP_sigma1 = phat_rl2(1);
AP_ccdf_rl = raylcdf(X_AP,AP_sigma1,'upper');

AP_lnmu = phat_ln2(1);
AP_lnsig = phat_ln2(2);
AP_ccdf_ln = logncdf(X_AP,AP_lnmu,AP_lnsig,'upper');

%==Extract GoF Data==%
AP_z_ex = expcdf(test_data2,AP_lambda1);
AP_z_gm = gamcdf(test_data2,AP_a1,AP_b1);
AP_z_rl = raylcdf(test_data2,AP_sigma1);
AP_z_ln = logncdf(test_data2,AP_lnmu,AP_lnsig);

AP_stats_ex = testStatistics(test_data2,AP_z_ex);
AP_stats_gm = testStatistics(test_data2,AP_z_gm);
AP_stats_rl = testStatistics(test_data2,AP_z_rl);
AP_stats_ln = testStatistics(test_data2,AP_z_ln);

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
axis([-inf,inf,1E-2,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
hold off
apfilename = [dir_ref,'/ActivityPotential_MLE.png'];
print(apfilename,'-dpng')
close(AP_fig);

%==Build and Return Relevant Data==%
times_struc_ex = struct('Scale',times_lambda);
times_struc_wb = struct('Scale',times_a,'Shape',times_b);

ap_struc_ex = struct('Scale',AP_lambda1);
ap_struc_gm = struct('Shape',AP_a1,'Scale',AP_b1);
ap_struc_rl = struct('Scale',AP_sigma1);
ap_struc_ln = struct('Location',AP_lnmu,'Scale',AP_lnsig);

times_size = size(times,2);
ap_size = size(activityPot,2);

times_pvals_ex = pvals_ex(times_size,times_lambda,times_stats_ex,3,6);
times_pvals_wb = pvals_wb(times_size,times_a,times_b,times_stats_wb,3,6);

ap_pvals_ex = pvals_ex(ap_size,AP_lambda1,AP_stats_ex,0,6);
ap_pvals_gm = pvals_gm(ap_size,AP_a1,AP_b1,AP_stats_gm,0,6);
ap_pvals_rl = pvals_rl(ap_size,AP_sigma1,AP_stats_rl,0,6);
ap_pvals_ln = pvals_ln(ap_size,AP_lnmu,AP_lnsig,AP_stats_ln,0,6);

times_EX = struct('Parameters',times_struc_ex,'Statistics',times_stats_ex,'pValues',times_pvals_ex);
times_WB = struct('Parameters',times_struc_wb,'Statistics',times_stats_wb,'pValues',times_pvals_wb);

AP_EX = struct('Parameters',ap_struc_ex,'Statistics',AP_stats_ex,'pValues',ap_pvals_ex);
AP_GM = struct('Parameters',ap_struc_gm,'Statistics',AP_stats_gm,'pValues',ap_pvals_gm);
AP_RL = struct('Parameters',ap_struc_rl,'Statistics',AP_stats_rl,'pValues',ap_pvals_rl);
AP_LN = struct('Parameters',ap_struc_ln,'Statistics',AP_stats_ln,'pValues',ap_pvals_ln);

IT_Struc = struct('Exponential',times_EX,'Weibull',times_WB);
AP_Struc = struct('Exponential',AP_EX,'Gamma',AP_GM,'Rayleigh',AP_RL,'LogNormal',AP_LN);
end