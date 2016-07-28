function [insert] = clineLaTex(paracount,height)

insert = cell(height,1);

if height == 1
    insert{1} = '\hline';
elseif height == 2
    text = cell(1,length(parcount));
    parfor i=1:length(parcount)
        if paracount(i) > 1
            text{i} = ['\cline{',2i+1,',',2i+2,'}'];
        end
    end
    insert{1} = [text{:}];
    insert{2} = '\hline';
elseif height == 3
    text1 = cell(1,length(parcount));
    text2 = cell(1,length(parcount));
    parfor i=1:length(parcount)
        if paracount(i) > 1
            text1{i} = ['\cline{',2i+1,',',2i+2,'}'];
        end
        if paracount(i) > 1
            text2{i} = ['\cline{',2i+1,',',2i+2,'}'];
        end
    end
    insert{1} = [text1{:}];
    insert{2} = [text2{:}];
    insert{3} = '\hline';
end