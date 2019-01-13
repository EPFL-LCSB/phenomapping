function [model, rules2check] = prepGPR(model)

grRules = model.grRules;
rules2check = {};
prob = 0;

% preprocessing
for i = 1:length(model.rxns)
    if length(grRules{i}) > 1
        % check if GPR has OR or AND rule together (origin of conflicts)
        ORgrtmp = ~isempty(find(~cellfun(@isempty,regexpi(grRules{i},'or','match'))));
        ANDgrtmp = ~isempty(find(~cellfun(@isempty,regexpi(grRules{i},'and','match'))));
        if (ORgrtmp && ANDgrtmp)
            % check if GPR has parenthesis near AND (origin of conflicts)
            flagANDPF = isempty(strfind(grRules{i},strcat('and',{' '},'(')));
            flagANDPB = isempty(strfind(grRules{i},strcat(')',{' '},'and')));
            if (flagANDPF==0) || (flagANDPB==0)
                prob = prob + 1;
                fprintf('Changing GPR:\t%d\n',prob);
                % split GPR based on parenthesis
                tmpGPR = strrep(grRules{i},' ','');
                tmpGPR = strrep(tmpGPR,'(',' ');
                tmpGPR = strrep(tmpGPR,')',' ');
                tmpGenes = splitString(tmpGPR,' ');
                
                ORgrtmp2 = zeros(length(tmpGenes),1);
                ANDgrtmp2 = zeros(length(tmpGenes),1);
                tmpGenes2 = cell(length(tmpGenes),1);
                lengthtmpGenes2 = zeros(length(tmpGenes),1);
                connector = zeros(length(tmpGenes),1);
                
                % check each section of GPR for AND/OR and nb of genes
                for iti = 1:length(tmpGenes)
                    ORgrtmp2(iti) = ~isempty(find(~cellfun(@isempty,regexpi(tmpGenes{iti},'or','match'))));
                    ANDgrtmp2(iti) = ~isempty(find(~cellfun(@isempty,regexpi(tmpGenes{iti},'and','match'))));
                    tmp = strrep(tmpGenes{iti},'or',' ');
                    tmp = strrep(tmp,'and',' ');
                    if ~isequal(strrep(tmp,{' '},''),{''})
                        tmpGenes2{iti} = splitString(tmp,' ');
                    else
                        tmpGenes2{iti} = {''};
                    end
                    lengthtmpGenes2(iti) = length(tmpGenes2{iti});
                    connector(iti) = isempty(tmpGenes2{iti}{1,1});
                end
                if length(tmpGenes)>3.5
                    fprintf('Correct the following GPR manually!:\t%d\n',i);
                end
                
                newgrRule = {};
                iti = 1; % this loop will only correct a GPR with two sections
                if (ANDgrtmp2(iti+1)==1) && (connector(iti+1)==1)
                    for tt1 = 1:lengthtmpGenes2(iti)
                        for tt2 = 1:lengthtmpGenes2(iti+2)
                            tmpnewgrRule = strcat('(',tmpGenes2{iti}{tt1},{' '},...
                                'and',{' '},tmpGenes2{iti+2}{tt2},')');
                            if isempty(newgrRule)
                                newgrRule =  tmpnewgrRule;
                            elseif (ORgrtmp2(iti)==1 || ORgrtmp2(iti+2)==1)
                                newgrRule = strcat(newgrRule,{' '},'or',{' '},tmpnewgrRule);
                            end
                        end
                    end
                elseif (ANDgrtmp2(iti+1)==1) && (~connector(iti+1)==1) && (lengthtmpGenes2(iti+1)==1)
                    for tt1 = 1:lengthtmpGenes2(iti)
                        tmpnewgrRule = strcat('(',tmpGenes2{iti}{tt1},{' '},...
                            'and',{' '},tmpGenes2{iti+1}{1},')');
                        if isempty(newgrRule)
                            newgrRule =  tmpnewgrRule;
                        elseif (ORgrtmp2(iti)==1)
                            newgrRule = strcat(newgrRule,{' '},'or',{' '},tmpnewgrRule);
                        end
                    end
                end
                rules2check{prob,1} = grRules{i};
                grRules(i) =  newgrRule;
                rules2check{prob,2} = newgrRule;
            end
        end
    end
end
model.grRules = grRules;
end
