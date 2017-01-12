function [data2global,FitTool,MLE,Moments] = analyse_GlobalClusteringCoeff(data,dir_ref)

Ymin = 1E-2;

cutExtreme = 0;

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

clustering = zeros(1,num_times);

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
    adj2 = thisadj^2;
    adj3 = thisadj^3;
    adj2sum = sum(sum(adj2));
    contrip = adj2sum - trace(adj2);
    if contrip==0
        clustering(m) = 0;
    else
        clustering(m) = trace(adj3)/contrip;
    end
end

clustering_autocorr = clustering;

clustering(clustering==0) = [];

FitTool = buildStruc_ExpGamRayLN_FitTool(clustering,dir_ref,'GlobalClusteringCoeff','Global Clustering Coefficient',cutExtreme,Ymin);
MLE = buildStruc_ExpGamRayLN_MLE(clustering,dir_ref,'GlobalClusteringCoeff','Global Clustering Coefficient',cutExtreme,Ymin);
Moments = buildStruc_ExpGamRayLN_Moments(clustering,dir_ref,'GlobalClusteringCoeff','Global Clustering Coefficient',cutExtreme,Ymin);

data2global = clustering;

maxTime = (num_times-1)*contact_time;
T = linspace(0,maxTime,num_times);

hwin = 50;
Tmod = T;
Tmod((num_times+1-hwin):num_times) = [];
Tmod(1:hwin) = [];

MA = zeros(1,num_times-2*hwin);
parfor i=1:(num_times-2*hwin)
    upper = i+2*hwin;
    MA(i) = sum(clustering_autocorr(i:upper))/(2*hwin+1);
end

clusteringfig = figure();
hold on
plot(T,clustering_autocorr)
plot(Tmod,MA,'LineWidth',4)
xlabel('Time (s)');
ylabel('Global Clustering Coefficient');
hold off
imagefilename = [dir_ref,'/GlobalClusteringCoeff_Additional_MovingAverage.png'];
print(imagefilename,'-dpng')
close(clusteringfig);

autocorrfig = figure();
subplot(3,1,1);
autocorr(clustering_autocorr,length(clustering_autocorr)-1);
subplot(3,1,2);
smallerlag = min(round(length(clustering_autocorr)-1,-2)/4,length(clustering_autocorr)-1);
autocorr(clustering_autocorr,smallerlag);
subplot(3,1,3);
evensmallerlag = min(round(length(clustering_autocorr)-1,-2)/16,length(clustering_autocorr)-1);
autocorr(clustering_autocorr,evensmallerlag);
imagefilename = [dir_ref,'/GlobalClusteringCoeff_Additional_AutoCorrelation.png'];
print(imagefilename,'-dpng');
close(autocorrfig);

end