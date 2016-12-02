%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates master files from the template with different asymmetric
% or migration rates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% open the template
f = fopen('Asymmetry_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

% while there are still lines
while ~feof(f) 
    % read line by line
    temp_lines{end+1,1} = fgets(f);   
end

% close the template xml    
fclose(f);  


% define the range of migration and coalescent rate asymmetries
migrationBias = logspace(-2,0,1000);
coalescentBias = logspace(-2,0,1000);

% make the xmls with asymmetric migration rates but symmetric coalescent
% rates
for mig = 1 : length(migrationBias)
    filename = sprintf('asymmetry_mig_%.5f_master',migrationBias(mig));
    fname = sprintf('./master/%s.xml',filename);   % set the file name
    p = fopen(fname,'w');
    c_nr = 1;
    m_nr = 1;

    for l = 1 : length(temp_lines)
        if ~isempty(strfind(temp_lines{l},'insert_coalescent'));
            fprintf(p,'%s',strrep(temp_lines{l},'insert_coalescent','1.0'));
            c_nr = c_nr +  1;
        elseif  ~isempty(strfind(temp_lines{l},'insert_migration'));
            if m_nr==1
                migration_rates = migrationBias(mig)*2/(migrationBias(mig)+1)*1.0;
            else
                migration_rates = 2/(migrationBias(mig)+1)*1.0;
            end
            fprintf(p,'%s',strrep(temp_lines{l},'insert_migration',num2str(migration_rates)));
            m_nr = m_nr + 1;
        elseif  ~isempty(strfind(temp_lines{l},'insert_samples'));                
            for a = 1 : 100
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',10*rand);
                fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="0"/>\n');
                fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
            end
            for a = 1 : 100
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',10*rand);
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
% make the xmls with symmetric migration and asymmetric coalescent rates
for coal = 1 : length(coalescentBias)
    filename = sprintf('asymmetry_coal_%.5f_master',coalescentBias(coal));
    fname = sprintf('./master/%s.xml',filename);   % set the file name
    p = fopen(fname,'w');
    c_nr = 1;
    m_nr = 1;

    for l = 1 : length(temp_lines)
        if ~isempty(strfind(temp_lines{l},'insert_coalescent'));
            if c_nr==1
                coalescent_rates = coalescentBias(coal)*2/(coalescentBias(coal)+1)*1.0;
            else
                coalescent_rates = 2/(coalescentBias(coal)+1)*1.0;
            end
            fprintf(p,'%s',strrep(temp_lines{l},'insert_coalescent',num2str(coalescent_rates)));
            c_nr = c_nr +  1;
        elseif  ~isempty(strfind(temp_lines{l},'insert_migration'));
            fprintf(p,'%s',strrep(temp_lines{l},'insert_migration','1.0'));
            m_nr = m_nr + 1;
        elseif  ~isempty(strfind(temp_lines{l},'insert_samples'));                
            for a = 1 : 100
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',10*rand);
                fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="0"/>\n');
                fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
            end
            for a = 1 : 100
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f">\n',10*rand);
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
