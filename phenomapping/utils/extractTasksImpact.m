function allTasksImpact = extractTasksImpact(impactTasks, geneList)

if (nargin < 2)
    geneList = {'infesible model'};
end

allTasksImpact = cell(length(geneList),3);

for i = 1:size(impactTasks,1)
    allTasksImpact(i,1) = geneList(i);
    if ~isempty(impactTasks{i})
        id = impactTasks{i}(1,1);
        name = impactTasks{i}(1,2);
        for j = 2:size(impactTasks{i},1)
            id = strcat(id,'|', impactTasks{i}(j,1));
            name = strcat(name,'|', impactTasks{i}(j,2));
        end
        allTasksImpact(i,2) = id;
        allTasksImpact(i,3) = name;
    else
        allTasksImpact(i,2) = {''};
        allTasksImpact(i,3) = {''};
    end
end