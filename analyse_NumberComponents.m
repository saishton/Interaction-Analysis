function [FitTool,MLE,Moments] = analyse_NumberComponents(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

components = zeros(1,num_times);

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
[~,thisComp,~] = networkComponents(thisadj);
thisComp(thisComp==1)=[];
thisCompCount = length(thisComp);
components(m) = thisCompCount;
end

FitTool = buildStruc_ExpGamRayLN_FitTool(components,dir_ref,'NumberComponents','Number of Components');
MLE = buildStruc_ExpGamRayLN_MLE(components,dir_ref,'NumberComponents','Number of Components');
Moments = buildStruc_ExpGamRayLN_Moments(components,dir_ref,'NumberComponents','Number of Components');

end