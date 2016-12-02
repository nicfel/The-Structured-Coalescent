%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates master files from the template with different sampling biases as
% well as different symmetric migration rate ranging from 0.01 to 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the template
f = fopen('SamplingBias_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

% while there are still lines
while ~feof(f) 
    temp_lines{end+1,1} = fgets(f);   % read line
end
    
fclose(f);  % close the template xml

% define the sample numbers for deme 1
bias = [10 30 50];
% Set the different migration rates
migration = [1 0.1 0.01];

for mig = 1 : length(migration)
    for l = 1 : length(bias)
        loc0 = bias(l);
        for i = 1 : 100   % print reps number of files    
            filename = sprintf('samplingBias_%d_mig_%.2f_sam_%d_master',i,migration(mig),loc0);
            fname = sprintf('./master/%s.xml',filename);   % set the file name
            p = fopen(fname,'w');
            for l = 1 : length(temp_lines)
                if ~isempty(strfind(temp_lines{l},'insert_coalescent'));
                    fprintf(p,'%s',strrep(temp_lines{l},'insert_coalescent','1.0'));
                elseif  ~isempty(strfind(temp_lines{l},'insert_migration'));
                    migration_rate = migration(mig);
                    fprintf(p,'%s',strrep(temp_lines{l},'insert_migration',num2str(migration_rate)));
                elseif  ~isempty(strfind(temp_lines{l},'insert_samples'));

                    loc1 = 100-loc0;
                    for a = 1 : loc0
                        fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',25*rand);
                        fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="0"/>\n');
                        fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
                    end
                    for a = 1 : loc1
                        fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',25*rand);
                        fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="1"/>\n');
                        fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
                    end            
                elseif ~isempty(strfind(temp_lines{l},'insert_filename'));
                    fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',filename));
                else
                    fprintf(p,'%s',temp_lines{l});  % print line unchanged
                end
            end
            fclose(p); %close file again
        end
    end
end