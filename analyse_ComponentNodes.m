function [FitTool,MLE,Moments] = analyse_ComponentNodes(data,dir_ref)

num_times = size(unique(data(:,1)),1);
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
contact_time = 20;

rawCompSizes = zeros(num_times,num_people);

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
[~,thisCompSizes,~] = networkComponents(thisadj);
thisPadding = num_people - length(thisCompSizes);
thisCompSizes = [thisCompSizes zeros(1,thisPadding)];
rawCompSizes(m,:) = thisCompSizes;
end

compSizes = rawCompSizes(:)';
compSizes(compSizes==0) = [];

FitTool = buildStruc_ExpGamRayLN_FitTool(compSizes,dir_ref,'ComponentNodes','Fraction of Nodes per Component');
MLE = buildStruc_ExpGamRayLN_MLE(compSizes,dir_ref,'ComponentNodes','Fraction of Nodes per Component');
Moments = buildStruc_ExpGamRayLN_Moments(compSizes,dir_ref,'ComponentNodes','Fraction of Nodes per Component');

end