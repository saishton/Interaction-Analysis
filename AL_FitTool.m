function [Structure] = AL_FitTool(data,dir_ref)

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

%==Plot clustering data==%
maxTime = (num_times-1)*contact_time;
T = linspace(0,maxTime,num_times);

hwin = 50;
Tmod = T;
Tmod((num_times+1-hwin):num_times) = [];
Tmod(1:hwin) = [];

MA = zeros(1,num_times-2*hwin);
parfor i=1:(num_times-2*hwin)
    upper = i+2*hwin;
    MA(i) = sum(clustering(i:upper))/(2*hwin+1);
end

clusteringfig = figure();
hold on
plot(T,clustering)
plot(Tmod,MA,'LineWidth',4)
xlabel('Time (s)');
ylabel('Clustering Coefficient');
hold off
imagefilename = [dir_ref,'/GlobalClusteringCoeff.png'];
print(imagefilename,'-dpng')
close(clusteringfig);

%==Fit Distributions for Number of Active Links==%
[F_links,X_links] = ecdf(numlinks);
ccdf_links = 1-F_links;
Xrem = [X_links(1);X_links(end-2,end)];
X_links = X_links(2:end-3);
ccdf_links = ccdf_links(2:end-3);

links_fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
links_ft_ex = fittype('expcdf(x,lambda,''upper'')','options',links_fo_ex);
[links_cf_ex,links_gof_ex] = fit(X_links,ccdf_links,links_ft_ex);
links_cv_ex = coeffvalues(links_cf_ex);

links_fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[1 1]);
links_ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',links_fo_gm);
[links_cf_gm,links_gof_gm] = fit(X_links,ccdf_links,links_ft_gm);
links_cv_gm = coeffvalues(links_cf_gm);

links_fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
links_ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',links_fo_rl);
[links_cf_rl,links_gof_rl] = fit(X_links,ccdf_links,links_ft_rl);
links_cv_rl = coeffvalues(links_cf_rl);

links_fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[0 1]);
links_ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',links_fo_ln);
[links_cf_ln,links_gof_ln] = fit(X_links,ccdf_links,links_ft_ln);
links_cv_ln = coeffvalues(links_cf_ln);

%==Extract Parameters==%
links_lambda = links_cv_ex(1);
links_ccdf_ex = expcdf(X_links,links_lambda,'upper');

links_a = links_cv_gm(1);
links_b = links_cv_gm(2);
links_ccdf_gm = gamcdf(X_links,links_a,links_b,'upper');

links_sigma = links_cv_rl(1);
links_ccdf_rl = raylcdf(X_links,links_sigma,'upper');

links_lnmu = links_cv_ln(1);
links_lnsig = links_cv_ln(2);
links_ccdf_ln = logncdf(X_links,links_lnmu,links_lnsig,'upper');

%==Extract GoF Data==%
dataMod = numlinks(~ismember(numlinks,Xrem));
links_test_data = sort(dataMod)';

links_z_ex = expcdf(links_test_data,links_lambda);
links_z_gm = gamcdf(links_test_data,links_a,links_b);
links_z_rl = raylcdf(links_test_data,links_sigma);
links_z_ln = logncdf(links_test_data,links_lnmu,links_lnsig);

links_stats_ex = testStatistics(links_test_data,links_z_ex);
links_stats_gm = testStatistics(links_test_data,links_z_gm);
links_stats_rl = testStatistics(links_test_data,links_z_rl);
links_stats_ln = testStatistics(links_test_data,links_z_ln);

links_stats_ex.Root_MSE = links_gof_ex.rmse;
links_stats_gm.Root_MSE = links_gof_gm.rmse;
links_stats_rl.Root_MSE = links_gof_rl.rmse;
links_stats_ln.Root_MSE = links_gof_ln.rmse;
links_stats_ex.R_Squared = links_gof_ex.rsquare;
links_stats_gm.R_Squared = links_gof_gm.rsquare;
links_stats_rl.R_Squared = links_gof_rl.rsquare;
links_stats_ln.R_Squared = links_gof_ln.rsquare;

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
xlabel('Number of Active Links');
ylabel('CCDF');
axis([1E0,1E2,1E-4,1E0]);
legend('Data','Exponential','Gamma','Rayleigh','Log-Normal');
imagefilename = [dir_ref,'/LinkActivations_FitTool.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Build and Return Relevant Data==%
links_struc_ex = struct('Scale',links_lambda);
links_struc_gm = struct('Shape',links_a,'Scale',links_b);
links_struc_rl = struct('Scale',links_sigma);
links_struc_ln = struct('Location',links_lnmu,'Scale',links_lnsig);

EX = struct('Parameters',links_struc_ex,'Statistics',links_stats_ex);
GM = struct('Parameters',links_struc_gm,'Statistics',links_stats_gm);
RL = struct('Parameters',links_struc_rl,'Statistics',links_stats_rl);
LN = struct('Parameters',links_struc_ln,'Statistics',links_stats_ln);

Structure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN);
end