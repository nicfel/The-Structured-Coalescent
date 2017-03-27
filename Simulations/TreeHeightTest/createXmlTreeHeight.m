%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the script creates the xml files used to sampled the tree height from the
% suggested structured coalescent distribution. It also creates an MASTER
% xml were the tree heights are inferred by simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define the number of lineages per deme
lineages = [3 2 1];
% define the pairwise coalescent rate in each deme
coalrate = [1 2 4];
% define all the migration rates between demes
migrate = [1 2 0.1 0.3 1 1]*0.01;
%define all the sampling time of the lineages
samplingtimes = [0 1 2 1 0 1];

%Creates all the xmls needed to compare tree height estimates
% %% make basta files
f = fopen('treeheight_basta.xml','r');    % open the template

temp_lines = cell(0,0);  % variable to save the lines

while ~feof(f)  % while there are still lines
    temp_lines{end+1,1} = fgets(f);   % read line
end
    
fclose(f);  % close the template xml

name = 'treeheight';
for j = 1 : length(lineages)
    name = sprintf('%s_%d',name,lineages(j));
end
name = sprintf('%s_coal',name);
for j = 1 : length(coalrate)
    name = sprintf('%s_%.3f',name,coalrate(j));
end
name = sprintf('%s_mig',name);
for j = 1 : length(migrate)
    name = sprintf('%s_%.3f',name,migrate(j));
end


fname = sprintf('./xmls/%s_basta.xml',name);   % set the file name
p = fopen(fname,'w');
for l = 1 : length(temp_lines)
    if ~isempty(strfind(temp_lines{l},'insert_leafs'));
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                fprintf(p, '\t\t<sequence id="inv%d_%d" taxon="inv%d_%d" totalcount="4" value="??"/>\n',b-1,a-1,b-1,a-1);
            end
        end  
    elseif ~isempty(strfind(temp_lines{l},'insert_trait'));
        traits = '';
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                traits = sprintf('%s,inv%d_%d=%d',traits,b-1,a-1,a-1);
            end
        end
        traits(1)=[];
        fprintf(p,strrep(temp_lines{l},'insert_trait',traits));
    elseif ~isempty(strfind(temp_lines{l},'insert_times'));
        times = '';
        c=1;
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                times = sprintf('%s,inv%d_%d=%.0f',times,b-1,a-1,samplingtimes(c));
                c = c+1;
            end
        end
        times(1)=[];
        fprintf(p,strrep(temp_lines{l},'insert_times',times));
    elseif ~isempty(strfind(temp_lines{l},'insert_coalescentRate'));
        line = strrep(temp_lines{l},'insert_coalescentRate',strtrim(sprintf('%.2f ',1./coalrate)));
        fprintf(p,'%s',strrep(line,'insert_dimension',num2str(length(lineages))));
    elseif  ~isempty(strfind(temp_lines{l},'insert_migrationRate'));
        line = strrep(temp_lines{l},'insert_migrationRate',strtrim(sprintf('%.4f ',migrate)));
        fprintf(p,'%s',strrep(line,'insert_dimension',num2str(length(lineages)*(length(lineages)-1))));
    elseif ~isempty(strfind(temp_lines{l},'insert_filename'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',[name '_basta']));
    else
        fprintf(p,'%s',temp_lines{l});  % print line unchanged
    end
end
fclose(p); %close file again

%% make esco files, the xml files for the exact structured coalescent

f = fopen('treeheight_esco.xml','r');    % open the template

temp_lines = cell(0,0);  % variable to save the lines

while ~feof(f)  % while there are still lines
    temp_lines{end+1,1} = fgets(f);   % read line
end
    
fclose(f);  % close the template xml

fname = sprintf('./xmls/%s_esco.xml',name);   % set the file name
p = fopen(fname,'w');
for l = 1 : length(temp_lines)
    if ~isempty(strfind(temp_lines{l},'insert_leafs'));
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                fprintf(p, '\t\t<sequence id="inv%d_%d" taxon="inv%d_%d" totalcount="4" value="??"/>\n',b-1,a-1,b-1,a-1);
            end
        end  
    elseif ~isempty(strfind(temp_lines{l},'insert_times'));
        times = '';
        c=1;
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                times = sprintf('%s,inv%d_%d=%.0f',times,b-1,a-1,samplingtimes(c));
                c = c+1;
            end
        end
        times(1)=[];
        fprintf(p,strrep(temp_lines{l},'insert_times',times));
    elseif ~isempty(strfind(temp_lines{l},'insert_coalescentRate'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_coalescentRate',strtrim(sprintf('%.1f ',coalrate))));
    elseif  ~isempty(strfind(temp_lines{l},'insert_migrationRate'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_migrationRate',strtrim(sprintf('%.4f ',migrate))));      
    elseif ~isempty(strfind(temp_lines{l},'insert_dimension'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_dimension',num2str(length(lineages))));
    elseif ~isempty(strfind(temp_lines{l},'insert_filename'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',[name '_esco']));
    else
        fprintf(p,'%s',temp_lines{l});  % print line unchanged
    end
end
fclose(p); %close file again


%% convert esco to masco files



f = fopen(sprintf('xmls/%s_esco.xml',name),'r');
g = fopen(sprintf('xmls/%s_masco.xml',name),'w');
while ~feof(f)
    line = fgets(f);
    if ~isempty(strfind(line,'ExactStructuredCoalescent'))
        fprintf(g, '%s', strrep(line,'ExactStructuredCoalescent','Masco'));
    elseif ~isempty(strfind(line,'_esco.'))
        fprintf(g, '%s', strrep(line,'_esco.','_masco.'));
    elseif ~isempty(strfind(line,'exactDensity'))
        fprintf(g, '%s', strrep(line,'exactDensity','mascoDensity'));
    else
        fprintf(g, '%s', line);
    end
end
fclose(f);
fclose(g);

%% convert esco to lisco files, the structured coalescent assuming lineages
% to be independent



f = fopen(sprintf('xmls/%s_esco.xml',name),'r');
g = fopen(sprintf('xmls/%s_lisco.xml',name),'w');
while ~feof(f)
    line = fgets(f);
    if ~isempty(strfind(line,'ExactStructuredCoalescent'))
        fprintf(g, '%s', strrep(line,'ExactStructuredCoalescent','IndependentStructuredCoalescent'));
    elseif ~isempty(strfind(line,'_esco.'))
        fprintf(g, '%s', strrep(line,'_esco.','_lisco.'));
    elseif ~isempty(strfind(line,'exactDensity'))
        fprintf(g, '%s', strrep(line,'exactDensity','independentDensity'));
    else
        fprintf(g, '%s', line);
    end
end
fclose(f);
fclose(g);

%% convert esco to sisco files, the structured coalescent according to the 
% volz approximation

f = fopen(sprintf('xmls/%s_esco.xml',name),'r');
g = fopen(sprintf('xmls/%s_sisco.xml',name),'w');
while ~feof(f)
    line = fgets(f);
    if ~isempty(strfind(line,'ExactStructuredCoalescent'))
        fprintf(g, '%s', strrep(line,'ExactStructuredCoalescent"','ApproximateStructuredCoalescent"'));
    elseif ~isempty(strfind(line,'_esco.'))
        fprintf(g, '%s', strrep(line,'_esco.','_sisco.'));
    elseif ~isempty(strfind(line,'exactDensity'))
        fprintf(g, '%s', strrep(line,'exactDensity','approximateDensity'));
    else
        fprintf(g, '%s', line);
    end
end
fclose(f);
fclose(g);


%% make master file to directly simulate tree heights

f = fopen('treeheight_master.xml','r');    % open the template

temp_lines = cell(0,0);  % variable to save the lines

while ~feof(f)  % while there are still lines
    temp_lines{end+1,1} = fgets(f);   % read line
end
    
fclose(f);  % close the template xml

fname = sprintf('./xmls/%s_master.xml',name);   % set the file name
p = fopen(fname,'w');
c=1;
m=1;
for l = 1 : length(temp_lines)
    if ~isempty(strfind(temp_lines{l},'insert_lineage'));
        s=1;
        for a = 1 : length(lineages)
            for b = 1 : lineages(a)
                fprintf(p, '\t\t\t<lineageSeedMultiple spec=''MultipleIndividuals'' copies="1" time="%.0f">\n',samplingtimes(s));
                fprintf(p, '\t\t\t\t<population spec=''Population'' type=''@L'' location="%d"/>\n',a-1);
                fprintf(p, '\t\t\t</lineageSeedMultiple>\n');
                s=s+1;
            end
        end  
    elseif ~isempty(strfind(temp_lines{l},'insert_coalescence'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_coalescence',num2str(coalrate(c)/2)));
        c = c + 1;
    elseif  ~isempty(strfind(temp_lines{l},'insert_migration'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_migration',num2str(migrate(m))));
        m = m + 1;
    elseif ~isempty(strfind(temp_lines{l},'insert_dimension'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_dimension',num2str(length(lineages))));
    elseif ~isempty(strfind(temp_lines{l},'insert_filename'));
        fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',[name '_master']));
    else
        fprintf(p,'%s',temp_lines{l});  % print line unchanged
    end
end
fclose(p); %close file again

