%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script creates master files with different migration rates. The
% migration rates are symmetric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open the template xml
f = fopen('migrationRates_master.xml','r');

% make a copy of all lines in the template
temp_lines = cell(0,0);
while ~feof(f)  
    temp_lines{end+1,1} = fgets(f);   
end    
fclose(f); 

% Use 1000 different symmetric migration rates between 10^-5 and 10^-0
migration_rates = logspace(-5,0,1000);


for i = 1 : length(migration_rates);   
    filename = sprintf('migrationRates_%.8f_master',migration_rates(i));        
    fname = sprintf('master/%s.xml',filename);   % set the file name
    p = fopen(fname,'w');
    for l = 1 : length(temp_lines)
        if  ~isempty(strfind(temp_lines{l},'insert_migration'));
            fprintf(p,'%s',strrep(temp_lines{l},'insert_migration',sprintf('%.8f',migration_rates(i))));        
        elseif ~isempty(strfind(temp_lines{l},'fileName'));
            l_tmp1 = strsplit(temp_lines{l},'fileName');    % split the line
            l_tmp2 = strsplit(l_tmp1{2},'.');    % split the line
            l_tmp3 = l_tmp2{1}(3:end);    % split the line
            line = strrep(temp_lines{l},l_tmp3,...
                sprintf('%s',filename));
            fprintf(p,'%s',line);
        else
            fprintf(p,'%s',temp_lines{l});  % print line unchanged
        end
    end
    fclose(p); %close file again
end

