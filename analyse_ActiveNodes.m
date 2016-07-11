function [FitTool,MLE,Moments] = analyse_ActiveNodes(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

nodes = zeros(1,num_times);

parfor m=1:num_times
    thisactive = zeros(1,num_people);
    current_time = (m-1)*contact_time;
    for i=1:data_length
        test_time = data(i,1);
        if test_time==current_time
            person1 = data(i,2);
            person2 = data(i,3);
            thisactive(person1) = 1;
            thisactive(person2) = 1;
        end
    end
    nodes(m) = sum(thisactive)/num_people;
end

FitTool = buildStruc_ExpGamRayLN_FitTool(nodes,dir_ref,'ActiveNodes','Percentage of Nodes Active');
MLE = buildStruc_ExpGamRayLN_MLE(nodes,dir_ref,'ActiveNodes','Percentage of Nodes Active');
Moments = buildStruc_ExpGamRayLN_Moments(nodes,dir_ref,'ActiveNodes','Percentage of Nodes Active');

end