function [data2global,FitTool,MLE,Moments] = analyse_ActivityPotential(data,dir_ref)

Ymin = 1E-1;

cutExtreme = 3;

number_rows = size(data,1);
number_people = max([data(:,2); data(:,3)]);
contact_time = 20;

%==Sort Data by ID==%
[~, order] = sort(data(:,3));
partsorteddata = data(order,:);
[~, order] = sort(partsorteddata(:,2));
sorteddata = partsorteddata(order,:);

%==Find Interaction Times==%
j = 1;
interactions = zeros(1,number_people);
step_vector = [contact_time 0 0];
while j<number_rows+1
    ID1 = sorteddata(j,2);
    ID2 = sorteddata(j,3);
    interactions(ID1) = interactions(ID1)+1;
    interactions(ID2) = interactions(ID2)+1;
    current_row = sorteddata(j,:);
    if j == number_rows
        next_row = [0 0 0]; 
    else
        next_row = sorteddata(j+1,:);
    end
    while isequal(next_row,current_row+step_vector)
        j = j+1;
        current_row = sorteddata(j,:);
        if j == number_rows
            next_row = [0 0 0];
        else
            next_row = sorteddata(j+1,:);
        end
    end
    j = j+1;
end
activityPot = 2*interactions/sum(interactions);

FitTool = buildStruc_ExpGamRayLN_FitTool(activityPot,dir_ref,'ActivityPotential','Activity Potential',cutExtreme,Ymin);
MLE = buildStruc_ExpGamRayLN_MLE(activityPot,dir_ref,'ActivityPotential','Activity Potential',cutExtreme,Ymin);
Moments = buildStruc_ExpGamRayLN_Moments(activityPot,dir_ref,'ActivityPotential','Activity Potential',cutExtreme,Ymin);

data2global = activityPot;
end