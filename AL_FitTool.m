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
Xrem = [X_links(end-2:end)];
X_links = X_links(2:end-3);
ccdf_links = ccdf_links(2:end-3);
dataMod = numlinks(~ismember(numlinks,Xrem));
links_test_data = sort(dataMod)';
maxX = max(X_links);

links_fo_po = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1]);
links_ft_po = fittype('poisscdf(x,lambda,''upper'')','options',links_fo_po);
[links_cf_po,links_gof_po] = fit(X_links,ccdf_links,links_ft_po);
links_cv_po = coeffvalues(links_cf_po);

%BI FT does not work

links_fo_nb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf 1],'StartPoint',[1 0.5]);
links_ft_nb = fittype('nbincdf(x,n,p,''upper'')','options',links_fo_nb);
[links_cf_nb,links_gof_nb] = fit(X_links,ccdf_links,links_ft_nb);
links_cv_nb = coeffvalues(links_cf_nb);

links_fo_ge = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[1],'StartPoint',[0.5]);
links_ft_ge = fittype('geocdf(x,r,''upper'')','options',links_fo_ge);
[links_cf_ge,links_gof_ge] = fit(X_links,ccdf_links,links_ft_ge);
links_cv_ge = coeffvalues(links_cf_ge);

%==Extract Parameters==%
links_pois_lambda = links_cv_po(1);
links_ccdf_po = poisscdf(X_links,links_pois_lambda,'upper');

%BI FT does not work

links_nb_p = links_cv_nb(2);
links_nb_r = links_cv_nb(1);
links_ccdf_nb = nbincdf(X_links,links_nb_r,links_nb_p,'upper');

links_geo_p = links_cv_ge(1);
links_ccdf_ge = geocdf(X_links,links_geo_p,'upper');

%==Extract GoF Data==%
links_z_po = poisscdf(links_test_data,links_pois_lambda);
links_z_nb = nbincdf(links_test_data,links_nb_r,links_nb_p);
links_z_ge = geocdf(links_test_data,links_geo_p);

links_stats_po = testStatistics(links_test_data,links_z_po);
links_stats_nb = testStatistics(links_test_data,links_z_nb);
links_stats_ge = testStatistics(links_test_data,links_z_ge);

links_stats_po.Root_MSE = links_gof_po.rmse;
links_stats_nb.Root_MSE = links_gof_nb.rmse;
links_stats_ge.Root_MSE = links_gof_ge.rmse;
links_stats_po.R_Squared = links_gof_po.rsquare;
links_stats_nb.R_Squared = links_gof_nb.rsquare;
links_stats_ge.R_Squared = links_gof_ge.rsquare;

%==Plotting==%
linksactivefig = figure();
hold on
plot(X_links,ccdf_links,'o');
plot(X_links,links_ccdf_po);
plot(X_links,links_ccdf_nb);
plot(X_links,links_ccdf_ge);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of Active Links');
ylabel('CCDF');
axis([-inf,inf,1E-4,1E0]);
legend('Data','Poisson','Neg. Binomial','Geometric');
imagefilename = [dir_ref,'/LinkActivations_Moments.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Build and Return Relevant Data==%
links_struc_po = struct('Mean',links_pois_lambda);
links_struc_nb = struct('Successes',links_nb_r,'Probability',links_nb_p);
links_struc_ge = struct('Probability',links_geo_p);
links_pvals_po = pvals_po(num_times,links_pois_lambda,links_stats_po,3,6);
links_pvals_nb = pvals_nb(num_times,links_nb_r,links_nb_p,links_stats_nb,3,6);
links_pvals_ge = pvals_ge(num_times,links_geo_p,links_stats_ge,3,6);

PO = struct('Parameters',links_struc_po,'Statistics',links_stats_po,'pValues',links_pvals_po);
NB = struct('Parameters',links_struc_nb,'Statistics',links_stats_nb,'pValues',links_pvals_nb);
GE = struct('Parameters',links_struc_ge,'Statistics',links_stats_ge,'pValues',links_pvals_ge);

Structure = struct('Poisson',PO,'NegBinomial',NB,'Geometric',GE);
end