function [data2global,FitTool,MLE,Moments] = analyse_ActiveEdges(data,dir_ref)

Ymin = 1E-3;

cutExtreme = 3;

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

links = zeros(1,num_times);
maxlinks = num_people*(num_people-1)/2;

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
    adjsum = sum(sum(thisadj));
    numlinks = adjsum/2;
    links(m) = numlinks/maxlinks;
end

FitTool = buildStruc_ExpGamRayLN_FitTool(links,dir_ref,'ActiveEdges','Fraction of Edges Active',cutExtreme,Ymin);
MLE = buildStruc_ExpGamRayLN_MLE(links,dir_ref,'ActiveEdges','Fraction of Edges Active',cutExtreme,Ymin);
Moments = buildStruc_ExpGamRayLN_Moments(links,dir_ref,'ActiveEdges','Fraction of Edges Active',cutExtreme,Ymin);

data2global = links;
end