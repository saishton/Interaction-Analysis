function [Structure] = NCT_MLE(data,dir_ref)

step = 20;
min_time = min(data(:,1));
max_time = max(data(:,1));
times = ((max_time-min_time)/step)+1;
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
rawactivity = zeros(data_length,num_people+1);

parfor i=1:data_length
    thisrawactivity = zeros(1,num_people+1);
    thisrawactivity(1) = data(i,1);
    person1 = data(i,2);
    person2 = data(i,3);
    thisrawactivity(person1+1) = 1;
    thisrawactivity(person2+1) = 1;
    rawactivity(i,:) = thisrawactivity;
end

activity = zeros(times,num_people);

for i=1:times
    currenttime = ((i-1)*step)+min_time;
    activerows = rawactivity(rawactivity(:,1)==currenttime,:);
    activerows = activerows(:,2:end);
    thisactivity = sum(activerows,1);
    thisactivity = (thisactivity>0);
    activity(i,:) = thisactivity;
end

activity = [activity; ones(1,num_people)];

long = activity(:);
long = long';
dlong = diff([1 long 1]);
startIndex = find(dlong < 0);
endIndex = find(dlong > 0)-1;
NoContact = endIndex-startIndex+1;
NoContact = NoContact*20;

%==Fit Distributions==%
[F_NCT,X_NCT] = ecdf(NoContact);
ccdf_NCT = 1-F_NCT;
Xrem = [X_NCT(1)];
X_NCT = X_NCT(2:end);
ccdf_NCT = ccdf_NCT(2:end);

dataMod = NoContact(~ismember(NoContact,Xrem));
test_data = sort(dataMod)';

mod_test_data = test_data;
mod_test_data(mod_test_data==0)=1E-99;

phat_ex = mle(mod_test_data,'distribution','Exponential');
phat_gm = mle(mod_test_data,'distribution','Gamma');
phat_rl = mle(mod_test_data,'distribution','Rayleigh');
phat_ln = mle(mod_test_data,'distribution','Lognormal');

%==Extract Parameters==%
NCT_lambda = phat_ex(1);
NCT_ccdf_ex = expcdf(X_NCT,NCT_lambda,'upper');

NCT_a = phat_gm(1);
NCT_b = phat_gm(2);
NCT_ccdf_gm = gamcdf(X_NCT,NCT_a,NCT_b,'upper');

NCT_sigma = phat_rl(1);
NCT_ccdf_rl = raylcdf(X_NCT,NCT_sigma,'upper');

NCT_lnmu = phat_ln(1);
NCT_lnsig = phat_ln(2);
NCT_ccdf_ln = logncdf(X_NCT,NCT_lnmu,NCT_lnsig,'upper');

%==Extract GoF Data==%
NCT_z_ex = expcdf(test_data,NCT_lambda);
NCT_z_gm = gamcdf(test_data,NCT_a,NCT_b);
NCT_z_rl = raylcdf(test_data,NCT_sigma);
NCT_z_ln = logncdf(test_data,NCT_lnmu,NCT_lnsig);

NCT_stats_ex = testStatistics(test_data,NCT_z_ex);
NCT_stats_gm = testStatistics(test_data,NCT_z_gm);
NCT_stats_rl = testStatistics(test_data,NCT_z_rl);
NCT_stats_ln = testStatistics(test_data,NCT_z_ln);

%==Plotting==%
NCT_fig = figure();
hold on
plot(X_NCT,ccdf_NCT,'o')
plot(X_NCT,NCT_ccdf_ex)
plot(X_NCT,NCT_ccdf_gm)
plot(X_NCT,NCT_ccdf_rl)
plot(X_NCT,NCT_ccdf_ln)
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Time Between Contacts (s)');
ylabel('CCDF');
axis([-inf,inf,1E-5,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Location','southwest');
hold off
imagefilename = [dir_ref,'/TimeBetweenContacts_MLE.png'];
print(imagefilename,'-dpng')
close(NCT_fig);

%==Build and Return Relevant Data==%
NCT_struc_ex = struct('Scale',NCT_lambda);
NCT_struc_gm = struct('Shape',NCT_a,'Scale',NCT_b);
NCT_struc_rl = struct('Scale',NCT_sigma);
NCT_struc_ln = struct('Location',NCT_lnmu,'Scale',NCT_lnsig);

NCT_size = size(NoContact,2);

NCT_pvals_ex = pvals_ex(NCT_size,NCT_lambda,NCT_stats_ex,0,3);
NCT_pvals_gm = pvals_gm(NCT_size,NCT_a,NCT_b,NCT_stats_gm,0,3);
NCT_pvals_rl = pvals_rl(NCT_size,NCT_sigma,NCT_stats_rl,0,3);
NCT_pvals_ln = pvals_ln(NCT_size,NCT_lnmu,NCT_lnsig,NCT_stats_ln,0,3);

NCT_EX = struct('Parameters',NCT_struc_ex,'Statistics',NCT_stats_ex,'pValues',NCT_pvals_ex);
NCT_GM = struct('Parameters',NCT_struc_gm,'Statistics',NCT_stats_gm,'pValues',NCT_pvals_gm);
NCT_RL = struct('Parameters',NCT_struc_rl,'Statistics',NCT_stats_rl,'pValues',NCT_pvals_rl);
NCT_LN = struct('Parameters',NCT_struc_ln,'Statistics',NCT_stats_ln,'pValues',NCT_pvals_ln);

Structure = struct('Exponential',NCT_EX,'Gamma',NCT_GM,'Rayleigh',NCT_RL,'LogNormal',NCT_LN);
end