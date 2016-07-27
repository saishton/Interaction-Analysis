function [data2global,FitTool,MLE,Moments] = analyse_TimeBetweenContacts(data,dir_ref)

Ymin = 1E-4;

cutExtreme = 3;

step = 20;
min_time = min(data(:,1));
max_time = max(data(:,1));
times = ((max_time-min_time)/step)+1;
data_length = size(data(:,1),1);
num_people = max([data(:,2); data(:,3)]);
rawactivity = zeros(data_length,num_people+1);

parfor i=1:data_length
    thisrawactivity = zeros(1,num_people+1);
    thisrawactivity(1) = data(i,1);
    person1 = data(i,2);
    person2 = data(i,3);
    thisrawactivity(person1+1) = 1;
    thisrawactivity(person2+1) = 1;
    rawactivity(i,:) = thisrawactivity;
end

activity = zeros(times,num_people);

parfor i=1:times
    currenttime = ((i-1)*step)+min_time;
    activerows = rawactivity(rawactivity(:,1)==currenttime,:);
    activerows = activerows(:,2:end);
    thisactivity = sum(activerows,1);
    thisactivity = (thisactivity>0);
    activity(i,:) = thisactivity;
end

activity = [activity; ones(1,num_people)];

long = activity(:);
long = long';
dlong = diff([1 long 1]);
startIndex = find(dlong < 0);
endIndex = find(dlong > 0)-1;
nocontact = endIndex-startIndex+1;
nocontact = nocontact*20;

FitTool = buildStruc_ExpGamRayLN_FitTool(nocontact,dir_ref,'TimeBetweenContacts','Length of Time Between Contacts',cutExtreme,Ymin);
MLE = buildStruc_ExpGamRayLN_MLE(nocontact,dir_ref,'TimeBetweenContacts','Length of Time Between Contacts',cutExtreme,Ymin);
Moments = buildStruc_ExpGamRayLN_Moments(nocontact,dir_ref,'TimeBetweenContacts','Length of Time Between Contacts',cutExtreme,Ymin);

data2global = nocontact;
end