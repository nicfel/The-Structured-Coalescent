%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');


for i = 1 : length(tree_files)

    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    
    namechar = strsplit(tree_files(i).name,'_');
    migrate = str2double(namechar{4});
    
    
    % coalescing
    tree_tmp2 = regexprep(t{1}(end-1),'&type="L",location="(\d)",reaction="Coalescence",time=(\d*).(\d*)','');

    %migrating
    tree_tmp1 = regexprep(tree_tmp2,'&type="L",location="(\d)",reaction="Migration",time=(\d*).(\d*)','');

    
    % make the MASTER tree compatible with BEAST2
    % sampling
    tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
    tip_locs = regexp(tree_tmp1,'[&type="L",location="(\d)",time=(\d*)\.(\d*)\]','match');
     
    for j = 1 : length(tip_locs{1})
        tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' tip_locs{1,1}{j}(22) 'kickout']);
        tree_tmp1 = strrep(tree_tmp1,'kickout','');
    end

    tree_tmp = regexprep(tree_tmp1,'(\d*)loc_','inv$1loc_');
    
    tree = strrep(tree_tmp{1},'[]','');
    if ~isempty(strfind(tree,']'))
        b = strfind(tree,']');
        c = tree((b-50):(b+50));
        disp(tree_files(i).name)
    end

    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    
    print_tree = tree;
    
    % make tripletts of all runs with different random initial values
    for tr = 1 : 3    
        % make the xmls for the structcoal
        flog = strrep(tree_files(i).name,'master.tree',sprintf('%dlisco',tr));
        fname = sprintf('xmls/%s.xml',flog);
        f = fopen(fname,'w');


        fprintf(g,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.structuredCoalescent.distribution:beast.structuredCoalescent.logger:beast.structuredCoalescent.operator:beast.structuredCoalescent.model:beast.structuredCoalescent.operators:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">\n');

        fprintf(g,'\t<data id="sequences" name="alignment">\n');
        for j = 1 : length(leafnames)
            fprintf(g,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="??"/>\n',leafnames{j},leafnames{j});
        end
        fprintf(g,'\t</data>\n');

        fprintf(g,'\t<run id="mcmc" spec="MCMC" chainLength="100000">\n');
        fprintf(g,'\t\t<state id="state" storeEvery="5000">\n');
        fprintf(g,'\t\t\t<tree id="tree" name="stateNode"/>\n');
        % draw the initial coalescent rates
        coalini = lognrnd(log(0.1), 0.4,2,1);
        while max(coalini)>=1.0
            coalini = lognrnd(log(0.1), 0.4,2,1);
        end
        fprintf(g,'\t\t\t<parameter id="coalRates" name="stateNode">%f %f</parameter>\n', coalini(1), coalini(2));
        % draw the initial migration rates
        migini = lognrnd(log(0.05),0.4,1,2);
        while max(migini)>=1.0
            migini = lognrnd(log(0.05),0.4,1,2);
        end
        fprintf(g,'\t\t\t<parameter id="migRates" name="stateNode">%f %f</parameter>\n', migini(1), migini(2));
        fprintf(g,'\t\t</state>\n');
        fprintf(g,'\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
        fprintf(g,'\t\t\tinitial="@tree" taxa="@sequences"\n');
        fprintf(g,'\t\t\tIsLabelledNewick="true"\n');
        fprintf(g,'\t\t\tnewick="%s"/>\n',print_tree);
        fprintf(g,' \t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@coalRates">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Exponential"  mean="2"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migRates">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Exponential"  mean="1"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t\t<distribution id="coalescent" spec="IndependentStructuredCoalescent" dim="2">\n');
        fprintf(g,'\t\t\t\t\t<structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@tree"/>\n');
        fprintf(g,'\t\t\t\t\t<coalescentRates idref="coalRates"/>\n');
        fprintf(g,'\t\t\t\t\t<migrationRates idref="migRates"/>\n');
        fprintf(g,'\t\t\t\t\t<parameter id="stepSize" name="timeStep">0.001</parameter>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t</distribution>\n');
        fprintf(g,'\t\t<operator id="CoalescentScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@coalRates" scaleAllIndependently="false" weight="1.0"/>\n');   
        fprintf(g,'\t\t<operator id="MigrationScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@migRates" scaleAllIndependently="false" weight="1.0"/>\n');   
        fprintf(g,'\t\t<logger id="tracelog" fileName="%s.log" logEvery="200" model="@posterior" sanitiseHeaders="true" sort="smart">\n',flog);
        fprintf(g,'\t\t\t<log idref="migRates"/>\n');
        fprintf(g,'\t\t\t<log idref="coalRates"/>\n');
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t\t<log idref="prior"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t\t<logger id="screenlog" logEvery="1000">\n');
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t</run>\n');
        fprintf(g,'</beast>\n');
        fclose(f);
    end
end

%% make SISCO files
% get all the LISCO files in the xml directory
iscoxmls = dir('xmls/*lisco.xml');


% convert each file such that it uses SISCO
for i = 1 : length(iscoxmls)
    f = fopen(sprintf('xmls/%s',iscoxmls(i).name),'r');
    g = fopen(sprintf('xmls/%s',strrep(iscoxmls(i).name,'lisco','sisco')),'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line,'IndependentStructuredCoalescent'))
            fprintf(g, '%s', strrep(line,'IndependentStructuredCoalescent"','ApproximateStructuredCoalescent"'));
        elseif ~isempty(strfind(line,'lisco.log'))
            fprintf(g, '%s', strrep(line,'lisco.log','sisco.log'));
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f);
    fclose(g);
end   
    
