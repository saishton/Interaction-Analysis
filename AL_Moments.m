function [Structure] = AL_Moments(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

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
    adjsum = sum(sum(thisadj));
    numlinks(m) = adjsum/2;
end

%==Cut Extreme Data and Estimate Moments==%
[F_links,X_links] = ecdf(numlinks);
ccdf_links = 1-F_links;
Xrem = [X_links(end-2:end)];
X_links = X_links(2:end-3);
ccdf_links = ccdf_links(2:end-3);
dataMod = numlinks(~ismember(numlinks,Xrem));
links_test_data = sort(dataMod)';

M1 = mean(dataMod);
M2 = mean(dataMod.^2);

%==Estimate Parameters for Distributions on Number of Active Links==%
links_pois_lambda = M1;
links_ccdf_po = poisscdf(X_links,links_pois_lambda,'upper');

links_bin_p = (1+M1)-(M2/M1);
links_bin_n = M1/links_bin_p;
links_ccdf_bi = binocdf(X_links,links_bin_n,links_bin_p,'upper');

links_nb_p = (M1/(M1^2-M2))-1;
links_nb_r = M1*(1-links_nb_p)/links_nb_p;
links_ccdf_nb = nbincdf(X_links,links_nb_r,links_nb_p,'upper');

links_geo_p = 1/(M1+1);
links_ccdf_ge = geocdf(X_links,links_geo_p,'upper');

%==Extract GoF Data==%
links_z_po = poisscdf(links_test_data,links_pois_lambda);
links_z_bi = binocdf(links_test_data,links_bin_n,links_bin_p);
links_z_nb = nbincdf(links_test_data,links_nb_r,links_nb_p);
links_z_ge = geocdf(links_test_data,links_geo_p);

links_stats_po = testStatistics(links_test_data,links_z_po);
links_stats_bi = testStatistics(links_test_data,links_z_bi);
links_stats_nb = testStatistics(links_test_data,links_z_nb);
links_stats_ge = testStatistics(links_test_data,links_z_ge);

links_stats_po.Root_MSE = sqrt(mean((links_z_po-links_test_data).^2));
links_stats_bi.Root_MSE = sqrt(mean((links_z_bi-links_test_data).^2));
links_stats_nb.Root_MSE = sqrt(mean((links_z_nb-links_test_data).^2));
links_stats_ge.Root_MSE = sqrt(mean((links_z_ge-links_test_data).^2));

TotSS = sum((links_test_data-mean(links_test_data)).^2);

links_stats_po.R_Squared = 1-sum((links_test_data-links_z_po).^2)/TotSS;
links_stats_bi.R_Squared = 1-sum((links_test_data-links_z_bi).^2)/TotSS;
links_stats_nb.R_Squared = 1-sum((links_test_data-links_z_nb).^2)/TotSS;
links_stats_ge.R_Squared = 1-sum((links_test_data-links_z_ge).^2)/TotSS;

%==Plotting==%
linksactivefig = figure();
hold on
plot(X_links,ccdf_links,'o');
plot(X_links,links_ccdf_po);
plot(X_links,links_ccdf_bi);
plot(X_links,links_ccdf_nb);
plot(X_links,links_ccdf_ge);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of Active Links');
ylabel('CCDF');
axis([-inf,inf,1E-4,1E0]);
legend('Data','Poisson','Binomial','Neg. Binomial','Geometric');
imagefilename = [dir_ref,'/LinkActivations_Moments.png'];
print(imagefilename,'-dpng')
close(linksactivefig);

%==Build and Return Relevant Data==%
links_struc_po = struct('Mean',links_pois_lambda);
links_struc_bi = struct('Trials',links_bin_n,'Probability',links_bin_p);
links_struc_nb = struct('Successes',links_nb_r,'Probability',links_nb_p);
links_struc_ge = struct('Probability',links_geo_p);
links_pvals_po = pvals_po(num_times,links_pois_lambda,links_stats_po,3,6);
links_pvals_bi = pvals_bi(num_times,links_bin_n,links_bin_p,links_stats_bi,3,6);
links_pvals_nb = pvals_nb(num_times,links_nb_r,links_nb_p,links_stats_nb,3,6);
links_pvals_ge = pvals_ge(num_times,links_geo_p,links_stats_ge,3,6);

PO = struct('Parameters',links_struc_po,'Statistics',links_stats_po,'pValues',links_pvals_po);
BI = struct('Parameters',links_struc_bi,'Statistics',links_stats_bi,'pValues',links_pvals_bi);
NB = struct('Parameters',links_struc_nb,'Statistics',links_stats_nb,'pValues',links_pvals_nb);
GE = struct('Parameters',links_struc_ge,'Statistics',links_stats_ge,'pValues',links_pvals_ge);

Structure = struct('Poisson',PO,'Binomial',BI,'NegBinomial',NB,'Geometric',GE);
end