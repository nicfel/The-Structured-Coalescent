%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the AIV.txt file to a fasta file with the requirements for it to
% run with LISCO and SISCO, i.e. add the _# at the end
% The AIV.txt file originates from the de Maio et al., 2015 paper. I
% just used their AIV_HA_BASTA_host_asymmetric.xml and changed it such that
% I can easily read in the data into matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% open the text file
f = fopen('AIV.txt','r');

% Counter for the Data structure. This is mainly needed in cases where not
% all lineages are used. Here all sequences were used. It would allow tough
% not all classes to be specified
i = 1;

location = cell(0,0);
species = cell(0,0);

% run while there are still lines left in the txt file, i.e. read line by
% line
while ~feof(f)
    line = fgets(f);
    % split the line after whitespaces
    tmp  = strsplit(line);
    
    % split the sequence names as they were defined in de Maio et al., 2015
    tmp2 = strsplit(tmp{1},'|');
    date = tmp2{2};
    name = tmp2{1};
    tmp3 = strsplit(tmp2{1},'/');
    tmp4 = strsplit(tmp3{1},'_');
    tmp5 = strsplit(tmp3{end},'_');
    location{end+1,1} = tmp5{3};
    
    % needed as sanity check to see if the text file was properly prepared.
    uniques{i} = tmp4{2};
    state = 'none';
    if ~isempty(strfind(lower(name),'alaska'))
        state = 'alaska_0';
        location{end+1,1} = state;
   elseif ~isempty(strfind(name,'Saskatchewan')) ||...
                    ~isempty(strfind(name,'NorthDakota'))||...
                    ~isempty(strfind(name,'BritishColumbia'))||...
                    ~isempty(strfind(name,'Alberta'))||...
                    ~isempty(strfind(name,'Washington'))
        state = 'northwest_1';
        location{end+1,1} = state;
    elseif ~isempty(strfind(name,'NewBrunswick')) ||...
                    ~isempty(strfind(name,'NovaScotia'))
        state = 'northeast_2';
        location{end+1,1} = state;
    elseif ~isempty(strfind(name,'California'))||...
                    ~isempty(strfind(name,'Mexico'))
                    state = 'southeast_3'; % should be south west. For the plotting use south west
        location{end+1,1} = state;
    elseif ~isempty(strfind(name,'Delaware'))||...
                    ~isempty(strfind(name,'NewJersey'))||...
                    ~isempty(strfind(name,'Maryland'))||...
                    ~isempty(strfind(name,'Pennsylvania'))||...
                    ~isempty(strfind(name,'NorthCarolina'))
        state = 'eastcoast_4';
        location{end+1,1} = state;
     elseif ~isempty(strfind(name,'Illinois'))||...
                    ~isempty(strfind(name,'Ohio'))||...
                    ~isempty(strfind(name,'Wisconsin'))
        state = 'northmideast_5';
        location{end+1,1} = state;  
     elseif ~isempty(strfind(name,'Mississippi'))||...
                    ~isempty(strfind(name,'Missouri'))||...
                    ~isempty(strfind(name,'Nebraska'))||...
                    ~isempty(strfind(name,'Texas'))
        state = 'center_6';
        location{end+1,1} = state; 
    else
        location{end+1,1} = tmp5{3};
    end
    
    if ~strcmp(state,'none')
        Data(i).Header = sprintf('%s|%s|%s|%s',name,tmp5{3},date,state);
        Data(i).Sequence = tmp{2};
        i = i + 1;
    end
end

% get the different locations
diff_locations = unique(location);
nr_sam = zeros(length(diff_locations),1);
for i = 1 : length(diff_locations)
    for j = 1 : length(location)
        if strcmp(location{j},diff_locations{i})
            nr_sam(i) = nr_sam(i) + 1;
        end
    end
end

% first delete the eventually existing fasta file. This needs to be done
% because matlab would otherwise just attach the new fasta data to the
% previously existing one
delete('AIV.fasta')
fastawrite('AIV.fasta',Data)

