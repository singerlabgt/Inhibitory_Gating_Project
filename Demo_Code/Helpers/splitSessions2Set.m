function [splitsets, splitanimals, uniqSess] = splitSessions2Set(index)
%this function splits rec days into 'learning set's where each set is a continuous
%set of rec days (usually 3 days) the animal was exposed to the same type
%of novel environment across days
%NJ 11/03/21
%NJ edit 01/13/22 to use regexp and flexibly pick sessions with novel 3-day exposure
%NJ edit 03/23/22 to add no-novelty paradigim animals to the split lists 

uniqSess = unique(index(:,[1:2,7]), 'rows'); %unique combo of animalID + recday

%only split sessions with regular expression of day 1-3 for novel paradigm, skip if not
patterns = '123';
numString = join(num2str(uniqSess(:,3))');
exp = regexp(numString, patterns);
exp = [exp, length(uniqSess) + 1]; %added this to deal with the last set
nonmatch = regexp(numString, patterns, 'split');
nonmatch = find(~cellfun(@isempty,nonmatch));

splitsets = cell(length(exp)-1, 1); splitanimals = zeros(length(exp)-1,1);
for ii = 1:length(exp)-1
    splitsets{ii} = uniqSess(exp(ii):exp(ii)+length(patterns)-1, :);
    splitanimals(ii) = unique( splitsets{ii}(:,1));
end

%find leftover sessions that did not have novelty exposure (eg. fam only)
%to add to the splitset list 
for ii = 1:length(nonmatch)
    leftoversess = uniqSess(exp(nonmatch(ii)-1)+length(patterns):exp(nonmatch(ii))-1, :); 
    leftoveranimal = unique(leftoversess(:,1)); 
    for an = 1:length(leftoveranimal)
        splitsets{end + 1} = leftoversess(leftoversess(:,1) == leftoveranimal(an), :);
        splitanimals(end + 1) = leftoveranimal(an);
    end
    
end