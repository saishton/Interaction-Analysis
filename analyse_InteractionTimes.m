function [data2global,FitTool,MLE,Moments] = analyse_InteractionTimes(data,dir_ref)

Ymin = 1E-3;

cutExtreme = 3;

number_rows = size(data,1);
contact_time = 20;

%==Sort Data by ID==%
[~, order] = sort(data(:,3));
partsorteddata = data(order,:);
[~, order] = sort(partsorteddata(:,2));
sorteddata = partsorteddata(order,:);

%==Find Interaction Times==%
times = zeros(1,number_rows);
j = 1;
times_k = 1;
step_vector = [contact_time 0 0];
while j<number_rows+1
    contact = contact_time;
    current_row = sorteddata(j,:);
    if j == number_rows
        next_row = [0 0 0]; 
    else
        next_row = sorteddata(j+1,:);
    end
    while isequal(next_row,current_row+step_vector)
        contact = contact+contact_time;
        j = j+1;
        current_row = sorteddata(j,:);
        if j == number_rows
            next_row = [0 0 0];
        else
            next_row = sorteddata(j+1,:);
        end
    end
    times(times_k) = contact;
    j = j+1;
    times_k = times_k+1;
end
times(times_k:end) = [];

FitTool = buildStruc_ExpMLGPWei_FitTool(times,dir_ref,'InteractionTimes','Length of Interaction',cutExtreme,Ymin);
MLE = buildStruc_ExpMLGPWei_MLE(times,dir_ref,'InteractionTimes','Length of Interaction',cutExtreme,Ymin);
Moments = buildStruc_ExpMLGPWei_Moments(times,dir_ref,'InteractionTimes','Length of Interaction',cutExtreme,Ymin);

data2global = times;
end