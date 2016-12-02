%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads trees and log files of mtt, esco, lisco and sisco and looks at
% the inferred root state probability when one migration rate is
% held constant and the other is varied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse the marginal state probabilities estimated using MultiTypeTree.
% The migration rate in one direction is only increase after n-iterations
% of state inference. 

% basis of the filenames
filename = 'rootStateProbabilities';

% the number of leafs used 
nrLeafs = 6;

% open the mtt *.trees and *.log files
f_log = fopen(sprintf('out/%s_mtt.log',filename),'r');

% read in the logs
t_log = textscan(f_log,'%s','CommentStyle','#');fclose(f_log);

% convert migration rates and root states to doubles
migrationRate = str2double(t_log{1,1}(13:5:end));
rootstate = str2double(t_log{1,1}(15:5:end));

unique_migration_rates = unique(migrationRate);

% initialize variables
MTT.nodeProbs= zeros(0,1);
MTT.migRate = zeros(0,1);

% get the mean root state probabilities
for i = 2 : length(unique_migration_rates)
    roots = rootstate(migrationRate==unique_migration_rates(i));
    MTT.nodeProbs(end+1,1) = mean(roots);
    MTT.migRate(end+1,1) = unique_migration_rates(i);
end

%% Analyse the root state probabilities calculated
% using the exact structured coalescent 


f_trees = fopen(sprintf('out/%s_esco.trees',filename),'r');
f_log = fopen(sprintf('out/%s_esco.log',filename),'r');

% read the trees and log files until theactual data is reached
c = true;
while c
    line_tree = fgets(f_trees);
    if ~isempty(strfind(line_tree,'tree STATE'))
        c = false;
    end
end
c = true;
while c
    line_log = fgets(f_log);
    if isempty(strfind(line_log,'#'))
        % The next line will be the first data line
        line_log = fgets(f_log);
        c = false;
    end
end

Esco.migRate = zeros(0,1);
Esco.nodeProbs = zeros(0,2*nrLeafs-1);

% Actually read in data
while isempty(strfind(line_tree,'End;'))
    trees_tmp = strsplit(line_tree);
    log_tmp = strsplit(line_log,'\t');
    
    % Save the migration rates
    Esco.migRate(end+1,1) = str2double(log_tmp{2});
    
    % Get the state probabilities
    tree = trees_tmp{end-1};
    node_probs_chars = strsplit(tree,'{');
    % The 4th (matlab counting) 7th, 9th and 10th
    % are internal nodes
    node_probs = zeros(1,2*nrLeafs-1);
    for i = 2 : length(node_probs_chars)
        node_probs_tmp = strsplit(node_probs_chars{i},',');
        node_probs(i-1) = str2double(node_probs_tmp{1});
    end    
    
    Esco.nodeProbs(end+1,:) = node_probs;
    
    line_tree = fgets(f_trees);
    line_log = fgets(f_log);
end


fclose(f_trees);
fclose(f_log);

%% Analyse the marginal state probabilities calculated
% using Lisco

f_trees = fopen(sprintf('out/%s_lisco.trees',filename),'r');
f_log = fopen(sprintf('out/%s_lisco.log',filename),'r');

% read the trees and log files until theactual data is reached
c = true;
while c
    line_tree = fgets(f_trees);
    if ~isempty(strfind(line_tree,'tree STATE'))
        c = false;
    end
end
c = true;
while c
    line_log = fgets(f_log);
    if isempty(strfind(line_log,'#'))
        % The next line will be the first data line
        line_log = fgets(f_log);
        c = false;
    end
end

Lisco.migRate = zeros(0,1);
Lisco.nodeProbs = zeros(0,2*nrLeafs-1);

% Actually read in data
while isempty(strfind(line_tree,'End;'))
    trees_tmp = strsplit(line_tree);
    log_tmp = strsplit(line_log,'\t');
    
    % Save the migration rates
    Lisco.migRate(end+1,1) = str2double(log_tmp{2});
    
    % Get the state probabilities
    tree = trees_tmp{end-1};
    node_probs_chars = strsplit(tree,'{');
    % The 4th (matlab counting) 7th, 9th and 10th
    % are internal nodes
    node_probs = zeros(1,2*nrLeafs-1);
    for i = 2 : length(node_probs_chars)
        node_probs_tmp = strsplit(node_probs_chars{i},',');
        node_probs(i-1) = str2double(node_probs_tmp{1});
    end    
    
    Lisco.nodeProbs(end+1,:) = node_probs;
    
    line_tree = fgets(f_trees);
    line_log = fgets(f_log);
end


fclose(f_trees);
fclose(f_log);


%% Analyse the marginal state probabilities calculated
% using sisco 

f_trees = fopen(sprintf('out/%s_sisco.trees',filename),'r');
f_log = fopen(sprintf('out/%s_sisco.log',filename),'r');

% read the trees and log files until theactual data is reached
c = true;
while c
    line_tree = fgets(f_trees);
    if ~isempty(strfind(line_tree,'tree STATE'))
        c = false;
    end
end
c = true;
while c
    line_log = fgets(f_log);
    if isempty(strfind(line_log,'#'))
        % The next line will be the first data line
        line_log = fgets(f_log);
        c = false;
    end
end

Sisco.migRate = zeros(0,1);
Sisco.nodeProbs = zeros(0,2*nrLeafs-1);

% Actually read in data
while isempty(strfind(line_tree,'End;'))
    trees_tmp = strsplit(line_tree);
    log_tmp = strsplit(line_log,'\t');
    
    % Save the migration rates
    Sisco.migRate(end+1,1) = str2double(log_tmp{2});
    
    % Get the state probabilities
    tree = trees_tmp{end-1};
    node_probs_chars = strsplit(tree,'{');
    node_probs = zeros(1,2*nrLeafs-1);
    for i = 2 : length(node_probs_chars)
        node_probs_tmp = strsplit(node_probs_chars{i},',');
        node_probs(i-1) = str2double(node_probs_tmp{1});
    end        
    Sisco.nodeProbs(end+1,:) = node_probs;
    
    line_tree = fgets(f_trees);
    line_log = fgets(f_log);
end
fclose(f_trees);
fclose(f_log);


% print the analysis to file such that one can plot it later in R
f = fopen('rootStateProbabilities.txt','w');
fprintf(f, 'migrationrate\tMTT\tesco\tLisco\tSisco\n');
for i = 1 : length(Esco.migRate)-1
    fprintf('%f\t%f\n',MTT.migRate(i),Esco.migRate(i))
    fprintf(f, '%f\t%f\t%f\t%f\t%f\n',Esco.migRate(i),1-MTT.nodeProbs(i),...
        Esco.nodeProbs(i,end),Lisco.nodeProbs(i,end),Sisco.nodeProbs(i,end));
end
fclose(f);





