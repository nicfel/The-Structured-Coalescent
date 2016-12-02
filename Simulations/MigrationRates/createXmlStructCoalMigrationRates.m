%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script creates ESCO, LISCO and SISCO files from master trees. The
% xmls are used to infer most likely symmetric migration rates. The
% inference is stopped once the maximum posterior is reached
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tree_files = dir('master/*.tree');

% needed???
% migrationratechar = cell(0,0);

for i = 1 : length(tree_files)    
    % find the true migration rate
    tmp = strsplit(tree_files(i).name,'_');
    migrationrate = str2double(tmp{2});

    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    
    % make the master tree compatible with BEAST2
    % sampling
    tip_locs = regexp(t{1}(end-1),'[&type="L",location="(\d)",time=(\d*).(\d*)\]','match');
    tree_tmp3 = t{1}{end-1};
    
    for j = 1 : length(tip_locs{1})
        tree_tmp3 = strrep(tree_tmp3,tip_locs{1}{j},['loc_' tip_locs{1,1}{j}(22) 'kickout']);
        tree_tmp3 = strrep(tree_tmp3,'kickout','');
    end
    
    tree_tmp4 = regexprep(tree_tmp3,'(\d*)loc_','inv$1loc_');
    
    % coalescing
    tree_tmp2 = regexprep(tree_tmp4,'&type="L",location="(\d)",reaction="Coalescence",time=(\d*).(\d*)','');

    % migrating
    tree_tmp1 = regexprep(tree_tmp2,'&type="L",location="(\d)",reaction="Migration",time=(\d*).(\d*)','');
    
    tree_tmp = regexprep(tree_tmp1,'E[-](\d)]',']');
    tree = strrep(tree_tmp,'[]','');
    if ~isempty(strfind(tree,']'))
        disp(tree)
    end

    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    
    print_tree = tree;       
    fname = regexprep(tree_files(i).name,'master.tree','esco.xml');
    flog = regexprep(fname,'.xml','');

    f = fopen(['xmls/' fname],'w');
    
    fprintf(f,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.structuredCoalescent.distribution:beast.structuredCoalescent.logger:beast.structuredCoalescent.operator:beast.structuredCoalescent.model:beast.structuredCoalescent.operators:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" version="2.0">\n');
    
   	fprintf(f,'\t<data id="sequences" name="alignment">\n');
    for j = 1 : length(leafnames)
        fprintf(f,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="??"/>\n',leafnames{j},leafnames{j});
    end
    fprintf(f,'\t</data>\n');
    
    fprintf(f,'\t<run id="mcmc" spec="MCMC" chainLength="100000">\n');
    fprintf(f,'\t\t<state id="state" storeEvery="5000">\n');
    fprintf(f,'\t\t\t<tree id="tree" name="stateNode"/>\n');
    fprintf(f,'\t\t\t<parameter id="coalRates" name="stateNode">2 2</parameter>\n');
    fprintf(f,'\t\t\t<parameter id="migRates" name="stateNode">%.12f %.12f</parameter>\n',migrationrate/100,migrationrate/100);
    fprintf(f,'\t\t</state>\n');
    fprintf(f,'\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
    fprintf(f,'\t\t\tinitial="@tree" taxa="@sequences"\n');
    fprintf(f,'\t\t\tIsLabelledNewick="true"\n');
    fprintf(f,'\t\t\tnewick="%s"/>\n',print_tree);
    fprintf(f,' \t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
    fprintf(f,'\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
    fprintf(f,'\t\t\t</distribution>\n');
    fprintf(f,'\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
    fprintf(f,'\t\t\t\t<distribution id="coalescent" spec="ExactStructuredCoalescent" dim="2">\n');
    fprintf(f,'\t\t\t\t\t<structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@tree"/>\n');
    fprintf(f,'\t\t\t\t\t<coalescentRates idref="coalRates"/>\n');
    fprintf(f,'\t\t\t\t\t<migrationRates idref="migRates"/>\n');
    fprintf(f,'\t\t\t\t\t<parameter id="stepSize" name="timeStep">0.005</parameter>\n');
    fprintf(f,'\t\t\t\t</distribution>\n');
    fprintf(f,'\t\t\t</distribution>\n');
    fprintf(f,'\t\t</distribution>\n');
    fprintf(f,'\t\t<operator id="CoalescentScaler" spec="ABCOperator" scaleFactor="1.02" scaleAll="true" parameter="@migRates" weight="1.0"/>\n');
    fprintf(f,'\t\t<logger id="tracelog" fileName="%s.log" logEvery="1" model="@posterior" sanitiseHeaders="true" sort="smart">\n',flog);
    fprintf(f,'\t\t\t<log idref="migRates"/>\n');
    fprintf(f,'\t\t\t<log idref="posterior"/>\n');
    fprintf(f,'\t\t</logger>\n');
    fprintf(f,'\t\t<logger id="screenlog" logEvery="10">\n');
    fprintf(f,'\t\t\t<log idref="posterior"/>\n');
    fprintf(f,'\t\t\t<log spec="EndCriterion" posterior="@coalescent" lowerBound=".01"/>\n');
    fprintf(f,'\t\t</logger>\n');
    fprintf(f,'\t</run>\n');
    fprintf(f,'</beast>\n');
    fclose(f);
end



%% creates LISCO files from ESCO files
escoxmls = dir('xmls/*esco.xml');

for i = 1 : length(escoxmls)
    f = fopen(sprintf('xmls/%s',escoxmls(i).name),'r');
    g = fopen(sprintf('xmls/%s',strrep(escoxmls(i).name,'esco','lisco')),'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line,'ExactStructuredCoalescent'))
            fprintf(g, '%s', strrep(line,'ExactStructuredCoalescent','IndependentStructuredCoalescent'));
        elseif ~isempty(strfind(line,'_esco.'))
            fprintf(g, '%s', strrep(line,'_esco.','_lisco.'));
        elseif ~isempty(strfind(line,'exactDensity'))
            fprintf(g, '%s', strrep(line,'exactDensity','independentDensity'));
        elseif ~isempty(strfind(line,'name="timeStep">0.005'))
            fprintf(g, '%s', strrep(line,'name="timeStep">0.005','name="timeStep">0.01'));            
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f);
    fclose(g);
end

%% make SISCO xml files from ESCO files

escoxmls = dir('xmls/*esco.xml');

for i = 1 : length(escoxmls)
    f = fopen(sprintf('xmls/%s',escoxmls(i).name),'r');
    g = fopen(sprintf('xmls/%s',strrep(escoxmls(i).name,'esco','sisco')),'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line,'ExactStructuredCoalescent'))
            fprintf(g, '%s', strrep(line,'ExactStructuredCoalescent"','ApproximateStructuredCoalescent"'));
        elseif ~isempty(strfind(line,'id="migRates"'))
            % since the migration rates are likely to be underestimated,
            % start at lower migration rates
             mrate = regexp(line, '(\d)\.(\d*)','match');
             newstartrate = str2double(mrate{1})*10^-4;
             line = strrep(line,mrate{1},sprintf('%.12f',newstartrate));
             fprintf(g, '%s', line);             
        elseif ~isempty(strfind(line,'_esco.'))
            fprintf(g, '%s', strrep(line,'_esco.','_sisco.'));
        elseif ~isempty(strfind(line,'exactDensity'))
            fprintf(g, '%s', strrep(line,'exactDensity','approximateDensity'));
        elseif ~isempty(strfind(line,'name="timeStep">0.005'))
            fprintf(g, '%s', strrep(line,'name="timeStep">0.005','name="timeStep">0.01'));            
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f);
    fclose(g);
end
