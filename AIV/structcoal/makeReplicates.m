%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make replicates with different initial migration and coalescent rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 

% read in the template xml
f = fopen('AIV.xml','r');

% save all the lines as a cell
lines = cell(0,0);
while ~feof(f)
    lines{end+1,1} = fgets(f);
end

% close the template again
fclose(f);

% make 3 replicates
for i = 1 : 3
    name = sprintf('AIV_%dlisco',i);
    g = fopen(['individual/' name '.xml'],'w');
    for l = 1 : length(lines)
        if ~isempty(strfind(lines{l,1},'insert_coal'))
            fprintf(g, strrep(lines{l,1}, 'insert_coal',...
                strtrim(sprintf('%f ', lognrnd(log(1),0.1,7,1)))));
        elseif ~isempty(strfind(lines{l,1},'insert_mig'))
            fprintf(g, strrep(lines{l,1}, 'insert_mig',...
                strtrim(sprintf('%f ', lognrnd(log(0.5),0.1,42,1)))));
        elseif ~isempty(strfind(lines{l,1},'spec="IndependentStructuredCoalescent"'))
            fprintf(g, strrep(lines{l,1}, 'spec="IndependentStructuredCoalescent"',...
                'spec="IndependentStructuredCoalescent" proportionalTimeStep="true"'));
        elseif ~isempty(strfind(lines{l,1},'insert_name'))
            fprintf(g, strrep(lines{l,1}, 'insert_name',name));
        else
            fprintf(g, lines{l,1});
        end
    end
    fclose(g);
end

for i = 1 : 3
    name = sprintf('AIV_%dsisco',i);
    g = fopen(['individual/' name '.xml'],'w');
    for l = 1 : length(lines)
        if ~isempty(strfind(lines{l,1},'insert_coal'))
            fprintf(g, strrep(lines{l,1}, 'insert_coal',...
                strtrim(sprintf('%f ', lognrnd(log(1),0.1,7,1)))));
        elseif ~isempty(strfind(lines{l,1},'insert_mig'))
            fprintf(g, strrep(lines{l,1}, 'insert_mig',...
                strtrim(sprintf('%f ', lognrnd(log(0.5),0.1,42,1)))));
        elseif ~isempty(strfind(lines{l,1},'name="timeStep">0.2'))
            fprintf(g, strrep(lines{l,1}, 'name="timeStep">0.2','name="timeStep">0.01'));
        elseif ~isempty(strfind(lines{l,1},'insert_name'))
            fprintf(g, strrep(lines{l,1}, 'insert_name',name));
        elseif ~isempty(strfind(lines{l,1},'independentDensity'))
            fprintf(g, strrep(lines{l,1}, 'independentDensity','approximateDensity'));
        elseif ~isempty(strfind(lines{l,1},'IndependentStructuredCoalescent'))
            fprintf(g, strrep(lines{l,1}, 'IndependentStructuredCoalescent','ApproximateStructuredCoalescent'));
        else
            fprintf(g, lines{l,1});
        end
    end
    fclose(g);
end

for i = 1 : 3
    name = sprintf('AIV_%dmasco',i);
    g = fopen(['individual/' name '.xml'],'w');
    for l = 1 : length(lines)
        if ~isempty(strfind(lines{l,1},'insert_coal'))
            fprintf(g, strrep(lines{l,1}, 'insert_coal',...
                strtrim(sprintf('%f ', lognrnd(log(1),0.1,7,1)))));
        elseif ~isempty(strfind(lines{l,1},'insert_mig'))
            fprintf(g, strrep(lines{l,1}, 'insert_mig',...
                strtrim(sprintf('%f ', lognrnd(log(0.5),0.1,42,1)))));
        elseif ~isempty(strfind(lines{l,1},'name="timeStep">0.2'))
            fprintf(g, strrep(lines{l,1}, 'name="timeStep">0.2','name="timeStep">0.1'));            
        elseif ~isempty(strfind(lines{l,1},'spec="IndependentStructuredCoalescent"'))
            fprintf(g, strrep(lines{l,1}, 'spec="IndependentStructuredCoalescent"',...
                'spec="Masco" proportionalTimeStep="true"'));
        elseif ~isempty(strfind(lines{l,1},'independentDensity'))
            fprintf(g, strrep(lines{l,1}, 'independentDensity','mascoDensity'));

        elseif ~isempty(strfind(lines{l,1},'insert_name'))
            fprintf(g, strrep(lines{l,1}, 'insert_name',name));
        else
            fprintf(g, lines{l,1});
        end
    end
    fclose(g);
end