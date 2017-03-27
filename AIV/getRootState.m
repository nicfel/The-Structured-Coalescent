%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads tree files of LISCO and SISCO and looks at
% the inferred root state probability 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% define the number of different states
states = 7;
% the folder where the states.trees files are
folder = 'structcoal/individual';
% the different states
states_names = {'alaska', 'northwest', 'northeast', 'southeast',...
                'eastcoast', 'northmideast', 'center'};

% read in tree files
f_masco = fopen([folder '/AIV_1masco.states.trees'], 'r');

% read in trees
rootState_masco1 = zeros(0,states);
while ~feof(f_masco)
    line = fgets(f_masco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        % Only take trees after a  burn in of 10%
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>300000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_masco1(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_masco);clear f_lisco line tmp1 tmp2 s2 s1

% read in tree files
rootState_masco2 = zeros(0,states);

f_masco = fopen([folder '/AIV_2masco.states.trees'], 'r');
while ~feof(f_masco)
    line = fgets(f_masco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>300000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_masco2(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_masco);clear f_lisco line tmp1 tmp2 s2 s1
% read in tree files
rootState_masco3 = zeros(0,states);
f_masco = fopen([folder '/AIV_3masco.states.trees'], 'r');
while ~feof(f_masco)
    line = fgets(f_masco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>300000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_masco3(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_masco);clear f_lisco line tmp1 tmp2 s2 s1




% Get the root state probabilities inferred using LISCO



% read in tree files
f_lisco = fopen([folder '/AIV_1lisco.states.trees'], 'r');

% read in trees
rootState_lisco1 = zeros(0,states);
while ~feof(f_lisco)
    line = fgets(f_lisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        % Only take trees after a  burn in of 10%
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_lisco1(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_lisco);clear f_lisco line tmp1 tmp2 s2 s1

% read in tree files
rootState_lisco2 = zeros(0,states);

f_lisco = fopen([folder '/AIV_2lisco.states.trees'], 'r');
while ~feof(f_lisco)
    line = fgets(f_lisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_lisco2(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_lisco);clear f_lisco line tmp1 tmp2 s2 s1
% read in tree files
rootState_lisco3 = zeros(0,states);
f_lisco = fopen([folder '/AIV_3lisco.states.trees'], 'r');
while ~feof(f_lisco)
    line = fgets(f_lisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_lisco3(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_lisco);clear f_lisco line tmp1 tmp2 s2 s1

% Get the root state probabilities inferred using SISCO

% read in tree files
f_sisco = fopen([folder '/AIV_1sisco.states.trees'], 'r');
rootState_sisco1 = zeros(0,states);
while ~feof(f_sisco)
    line = fgets(f_sisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_sisco1(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_sisco);clear f_sisco line tmp1 tmp2 s2 s1

% read in tree files
f_sisco = fopen([folder '/AIV_2sisco.states.trees'], 'r');
rootState_sisco2 = zeros(0,states);
while ~feof(f_sisco)
    line = fgets(f_sisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_sisco2(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_sisco);clear f_sisco line tmp1 tmp2 s2 s1

% read in tree files
f_sisco = fopen([folder '/AIV_3sisco.states.trees'], 'r');
rootState_sisco3 = zeros(0,states);
while ~feof(f_sisco)
    line = fgets(f_sisco);
    tmp1 = strsplit(line);
    if length(tmp1)>1
        tmp2 = strsplit(tmp1{2},'_');
        if ~isempty(strfind(line, 'tree STATE')) && str2double(tmp2{2})>2000000
            s2 = strsplit(line,'}');
            s1 = strsplit(s2{end-1},'{');
            rootState_sisco3(end+1,:) = str2double(strsplit(s1{2},','));
        end
    end
end
fclose(f_sisco);clear f_sisco line tmp1 tmp2 s2 s1

% calculate the mean root state probabilities for the combined runs
m = mean([rootState_masco1;rootState_masco2;rootState_masco3]);
l = mean([rootState_lisco1;rootState_lisco2;rootState_lisco3]);
s = mean([rootState_sisco1;rootState_sisco2;rootState_sisco3]);

% calculate the mean root state probabilities for each individual run
m1 = mean(rootState_masco1)';
m2 = mean(rootState_masco2)';
m3 = mean(rootState_masco3)';
l1 = mean(rootState_lisco1)';
l2 = mean(rootState_lisco2)';
l3 = mean(rootState_lisco3)';
s1 = mean(rootState_sisco1)';
s2 = mean(rootState_sisco2)';
s3 = mean(rootState_sisco3)';


% Combine the individual runs plus the combined mean values
m_plot = [m1,m2,m3,m'];
l_plot = [l1,l2,l3,l'];
s_plot = [s1,s2,s3,s'];
d_plot = [m',l',s'];


% Plot the mean root state probabilities of the individual chains to see
% if they are consistent between runs. A relatively high ESS value of the
% combined chains for each parameter should be enough to ensure
% convergence to the same optima, but it can act as a sanity check.
figure()
subplot(4,1,1)
bar(m_plot)
title('individual chains: MASCO')
set(gca,'XTickLabel',states_names);
legend('1^{st} chain','2^{nd} chain','3^{rd} chain','combined',...
    'Location','north','Orientation','horizontal');

subplot(4,1,2)
bar(l_plot)
title('individual chains: LISCO')
set(gca,'XTickLabel',states_names);
legend('1^{st} chain','2^{nd} chain','3^{rd} chain','combined',...
    'Location','north','Orientation','horizontal');
subplot(4,1,3)
bar(s_plot)
title('individual chains: SISCO')
set(gca,'XTickLabel',states_names);
legend('1^{st} chain','2^{nd} chain','3^{rd} chain','combined',...
    'Location','north','Orientation','horizontal');
subplot(4,1,4)
bar(d_plot(1:7,:))
title('combined chains: MASCO LISCO and SISCO')
set(gca,'XTickLabel',states_names);
legend('MASCO','LISCO','SISCO',...
    'Location','north','Orientation','horizontal');




% write the mean root state probabilities to a files such that they can be
% plotted later on in R
f = fopen('rootStates.txt','w');

states_names = {'alaska', 'northwest', 'northeast', 'southeast',...
                'eastcoast', 'northmideast', 'center'};
fprintf(f,'Method\t%s\n',strtrim(sprintf('%s\t',states_names{:})));
fprintf(f,'masco\t%s\n',strtrim(sprintf('%f\t',m)));
% fprintf(f,'lisco\t%s\n',strtrim(sprintf('%f\t',l)));
fprintf(f,'sisco\t%s\n',strtrim(sprintf('%f\t',s)));
fclose(f);clear f ans


