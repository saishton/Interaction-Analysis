function [Analysis] = analyse_interactions(input_filename,structure,timestamp)

iF = 'input';
oF = ['output_',timestamp];

clean_input = strrep(input_filename, '.', '');
dir_ref = [oF,'\',clean_input];
mkdir(dir_ref);

input = [iF,'/',input_filename];

fid = fopen(input);
rawdata = textscan(fid,structure,'Delimiter',',');
fclose(fid);

%==Extract and Clean Data==%
data = cell2mat(rawdata);
data(:,1) = data(:,1)-data(1,1);
lowestID = min(min(data(:,2)),min(data(:,3)));
data(:,2) = data(:,2)-lowestID+1;
data(:,3) = data(:,3)-lowestID+1;
number_rows = size(data,1);
parfor i=1:number_rows
    thisrow = data(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data(i,:) = thisrow;
    end
end
all_IDs = [data(:,2); data(:,3)];
all_active = unique(all_IDs);
num_people = size(all_active,1);
data2 = data(:,2);
data3 = data(:,3);
for i=1:num_people
    oldID = all_active(i);
    data2(data2==oldID) = -i;
    data3(data3==oldID) = -i;
end
data(:,2) = -data2;
data(:,3) = -data3;

%==Fit Distributions to Number of Active Links==%
ActiveLinks_FitTool = AL_FitTool(data,dir_ref);
ActiveLinks_Moments = AL_Moments(data,dir_ref);
ActiveLinks_MLE = AL_MLE(data,dir_ref);

%==Fit Distributions to Interaction Times & Activity Potential==%
[IT_FitTool,AP_FitTool] = IT_AP_FitTool(data,dir_ref);
[IT_Moments,AP_Moments] = IT_AP_Moments(data,dir_ref);
[IT_MLE,AP_MLE] = IT_AP_MLE(data,dir_ref);

InteractionTimes_FitTool = IT_FitTool;
ActivityPotential_FitTool = AP_FitTool;
InteractionTimes_Moments = IT_Moments;
ActivityPotential_Moments = AP_Moments;
InteractionTimes_MLE = IT_MLE;
ActivityPotential_MLE = AP_MLE;

%==Fit Distributions for Lengths of Time without a Specific Individual in an Interaction==%
NoContactTimes_FitTool = NCT_FitTool(data,dir_ref);
NoContactTimes_Moments = NCT_Moments(data,dir_ref);
NoContactTimes_MLE = NCT_MLE(data,dir_ref);

%==Fit Distributions for Number of Active Nodes==%
NodesActive_FitTool = ActiveNodes_FitTool(data,dir_ref);
NodesActive_Moments = ActiveNodes_Moments(data,dir_ref);
NodesActive_MLE = ActiveNodes_MLE(data,dir_ref);

datafilename = [dir_ref,'/Distributions.mat'];
save(datafilename,...
    'ActiveLinks_FitTool',...
    'InteractionTimes_FitTool',...
    'ActivityPotential_FitTool',...
    'NoContactTimes_FitTool',...
    'NodesActive_FitTool',...
    'ActiveLinks_Moments',...
    'InteractionTimes_Moments',...
    'ActivityPotential_Moments',...
    'NoContactTimes_Moments',...
    'NodesActive_Moments',...
    'ActiveLinks_MLE',...
    'InteractionTimes_MLE',...
    'ActivityPotential_MLE',...
    'NoContactTimes_MLE',...
    'NodesActive_MLE'...
)

Analysis = struct(  'ActiveLinks_FitTool',ActiveLinks_FitTool,...
                    'ActiveLinks_Moments',ActiveLinks_Moments,...
                    'ActiveLinks_MLE',ActiveLinks_MLE,...
                    'InteractionTimes_FitTool',InteractionTimes_FitTool,...
                    'InteractionTimes_Moments',InteractionTimes_Moments,...
                    'InteractionTimes_MLE',InteractionTimes_MLE,...
                    'ActivityPotential_FitTool',ActivityPotential_FitTool,...
                    'ActivityPotential_Moments',ActivityPotential_Moments,...
                    'ActivityPotential_MLE',ActivityPotential_MLE,...
                    'NoContactTimes_FitTool',NoContactTimes_FitTool,...
                    'NoContactTimes_Moments',NoContactTimes_Moments,...
                    'NoContactTimes_MLE',NoContactTimes_MLE,...
                    'NodesActive_FitTool',NodesActive_FitTool,...
                    'NodesActive_Moments',NodesActive_Moments,...
                    'NodesActive_MLE',NodesActive_MLE,...
                    'NumberPeople',num_people);

analysis2latex(Analysis,dir_ref);
end