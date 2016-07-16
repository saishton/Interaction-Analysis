function [FitTool,MLE,Moments] = analyse_ComponentEdges(data,dir_ref)

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
[~,~,thisCompGroups] = networkComponents(thisadj);
thisNumComps = length(thisCompGroups);
thisEdges = zeros(1,thisNumComps);
    for j=1:thisNumComps
        thisNodes = cell2mat(thisCompGroups(j));
        thisNodesSize = length(thisNodes);
        if thisNodesSize == 1
            thisEdges(j)=0;
        else
            thisMaxEdges = thisNodesSize*(thisNodesSize-1)/2;
            thisSubMat = thisadj(thisNodes,thisNodes);
            thisAdjSum = sum(sum(thisSubMat));
            thisNumEdges = thisAdjSum/2;
            thisFracEdges = thisNumEdges/thisMaxEdges;
            thisEdges(j) = thisFracEdges;
        end
    end
thisPadding = num_people - length(thisEdges);
thisEdges = [thisEdges zeros(1,thisPadding)];
rawFracEdges(m,:) = thisEdges;
end

compEdgeFracs = rawFracEdges(:)';
compEdgeFracs(compEdgeFracs==0) = [];

FitTool = buildStruc_ExpGamRayLN_FitTool(compEdgeFracs,dir_ref,'ComponentEdges','Fraction of Edges Active per Component');
MLE = buildStruc_ExpGamRayLN_MLE(compEdgeFracs,dir_ref,'ComponentEdges','Fraction of Edges Active per Component');
Moments = buildStruc_ExpGamRayLN_Moments(compEdgeFracs,dir_ref,'ComponentEdges','Fraction of Edges Active per Component');

end