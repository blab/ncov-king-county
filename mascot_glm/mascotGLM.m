clear;
%start_date = '2020-01-31';
%clade = {'subsampled','all', 'Alpha', 'Epsilon', 'Gamma'};
clade = {'subsampled', 'weighted', 'random'}; 
%clade = {'weighted'};
fastafiles = dir('data/*.fasta');
for file_names = 1 : length(fastafiles)
    s = fastaread(['data/' fastafiles(file_names).name]);
    seq_id=cell(0,0);
    for i = 1 : length(s)
        seq_id{i} = s(1).Header;
    end
    for v = 1: length(clade)
        f = fopen(['data/' strrep(fastafiles(file_names).name, 'fasta','tsv')]);
        dat = textscan(f, '%s\t%s\t%s\t%s\t%s%s\n','HeaderLines',1);fclose(f);
        seq_names = cell(0,0);
        for j = 1 : length(s)
            seq_names{j,1} = s(j).Header;
         end
        all_clusters = unique(dat{6});
    
        subsample_n = 3000;           
        
        if strcmp(clade{v}, 'subsampled')
            
            
            date_list = year(datetime(dat{1,4}))*100+ week(datetime(dat{1,4}));

            u_date = unique(date_list);
            tmp_names = dat{1};
            sub_names = cell(0,0);
            counter = 1;
            while counter <= subsample_n
                samp_date = randsample(u_date,1);
                names = tmp_names(ismember(date_list,samp_date));
                if ~isempty(names) 
                    rand_name = randsample(names, 1);
                    if ~isempty(rand_name)
                        sub_names{counter,1} = rand_name;
                        ind_name = find(ismember(tmp_names, sub_names{counter,1}));
                        tmp_names(ind_name) = [];
                        date_list(ind_name) = [];
                        counter = counter + 1 ;
                    end
                else 

                end
            end 

            flat_sub_names = [sub_names{:}];
            for a = length(dat{1}):-1:1
                if sum(ismember(dat{1}{a}, flat_sub_names)) == 0
                    dat{1}(a) = [];
                    dat{2}(a) = [];
                    dat{3}(a) = [];
                    dat{4}(a) = [];
                    dat{5}(a) = [];
                    dat{6}(a) = [];
                end 
            end
        end
        
        
        if strcmp(clade{v}, 'random')
            
            tmp_names = [dat{1}];
            sub_names = [];
            sub_names = randsample(tmp_names, subsample_n);
            flat_sub_names = sub_names;
            for a = length(dat{1}):-1:1
                if sum(ismember(dat{1}{a}, flat_sub_names)) == 0
                    dat{1}(a) = [];
                    dat{2}(a) = [];
                    dat{3}(a) = [];
                    dat{4}(a) = [];
                    dat{5}(a) = [];
                    dat{6}(a) = [];
                end 
            end
        end
        
        
        if strcmp(clade{v}, 'weighted')
            
            %u_raw_date = unique(datetime(dat{1,4}));
            T = readtable('~/Desktop/gitrepos/ncov-king-county/analysis/data-files/kc_hosp_weekly.csv');
            available_dates = T.Admission_Date;
            available_date_list = year(available_dates)*100+ week(available_dates);
            hosp = T.Hospitalizations;
            

            
            
            date_list = year(datetime(dat{1,4}))*100+ week(datetime(dat{1,4}));

            u_date = unique(date_list);
            weight_hosp = hosp(ismember(available_date_list, u_date));
            
            full_year_week = available_date_list(ismember(available_date_list, u_date));
            tmp_names = dat{1};
            sub_names = cell(0,0);
            counter = 1;
            while counter <= subsample_n
                samp_date = randsample(full_year_week,1, true, weight_hosp);
                names = tmp_names(ismember(date_list,samp_date));
                if ~isempty(names) 
                    rand_name = randsample(names, 1);
                    if ~isempty(rand_name)
                        sub_names{counter,1} = rand_name;
                        ind_name = find(ismember(tmp_names, sub_names{counter,1}));
                        tmp_names(ind_name) = [];
                        date_list(ind_name) = [];
                        counter = counter + 1 ;
                    end
                else 

                end
            end 

            flat_sub_names = [sub_names{:}];
            for a = length(dat{1}):-1:1
                if sum(ismember(dat{1}{a}, flat_sub_names)) == 0
                    dat{1}(a) = [];
                    dat{2}(a) = [];
                    dat{3}(a) = [];
                    dat{4}(a) = [];
                    dat{5}(a) = [];
                    dat{6}(a) = [];
                end 
            end
        end
        

        all_dates = dat{4};
        max_date = max(datenum(all_dates, 'yyyy-mm-dd'));
        min_date = min(datenum(all_dates, 'yyyy-mm-dd'));
        adjusted_start_date = datetime(min_date,'ConvertFrom','datenum') - calmonths(2);
        rate_shift = [1/366:1/366:(max_date-datenum(adjusted_start_date))/366];
        rate_shift_immi = [14/366:14/366:(max_date-datenum(adjusted_start_date))/366, 2.5];
        
        mcov_forward = repmat([1,0], 1, length(rate_shift)) ;
        mcov_back = repmat([0,1], 1, length(rate_shift)) ;

        local_clusters = unique(dat{6});
        writecell([dat{:}],['clusters' clade{v}  '.csv'])
        f = fopen(['template_glm_from_multicoal_nomig_' strrep(fastafiles(file_names).name, 'fasta','xml')]);
        g = fopen(['xmls/TESTglm_mcc_map_' clade{v} strrep(fastafiles(file_names).name, 'fasta','xml')],'w');
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_data')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t<data id="S%s" spec="Alignment" name="alignment">\n', local_clusters{lc});                
                    names = dat{1}(ismember(dat{6},local_clusters{lc}));
                        for k = 1 :length(names)
                            ind_name = find(ismember(seq_names, names{k}));
                            fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',...
                                s(ind_name).Header, s(ind_name).Header, s(ind_name).Sequence);
                        end
                        fprintf(g, '\t</data>\n');

                end

            elseif contains(line, 'insert_tree')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t<tree id="Tree.t:S%s" spec="beast.evolution.tree.Tree" name="stateNode">\n', local_clusters{lc});                            
                    dates = 'rem';
                    names = dat{1}(ismember(dat{6},local_clusters{lc}));
                    for k = 1 :length(names)
                        raw_date = dat{4}{ismember(dat{1},names{k})};
                        split_date = strsplit(raw_date,'-');
                        dates = [dates ',' names{k} '=' num2str(str2double(split_date{1}) +(datenum(raw_date, 'yyyy-mm-dd') - datenum(split_date{1}, 'yyyy'))/(datenum(num2str(str2double(split_date{1})+1), 'yyyy') - datenum(split_date{1}, 'yyyy'))) ];
                    end
                    dates = strrep(dates, 'rem,','');
                    fprintf(g, '\t\t\t<trait id="dateTrait.t:S%s" spec="beast.evolution.tree.TraitSet" dateFormat="decimal" traitname="date" value="%s">\n', local_clusters{lc}, dates);
                    fprintf(g, '\t\t\t\t<taxa id="TaxonSet.S%s" spec="TaxonSet">\n', local_clusters{lc});
                    fprintf(g, '\t\t\t\t\t<alignment idref="S%s"/>\n', local_clusters{lc});
                    fprintf(g, '\t\t\t\t</taxa>\n');
                    fprintf(g, '\t\t\t</trait>\n');
                    fprintf(g, '\t\t\t<taxonset idref="TaxonSet.S%s"/>\n', local_clusters{lc});
                    fprintf(g, '\t\t</tree>\n');
              %      fprintf(g, '\t\t<parameter id="rootLength:S%s" spec="parameter.RealParameter" lower="0.0" upper="0.25" name="stateNode">0.01</parameter>\n', local_clusters{lc});
                end

            elseif contains(line, 'insert_parameters')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t<parameter id="rootLength:S%s" name="stateNode" upper="1.0" dimension="1">0.001</parameter>\n', local_clusters{lc});
                end
                 
                 fprintf(g,'\t\t<parameter id="rateShifts.immi" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shift_immi));

            elseif contains(line, 'insert_init_tree')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t <init id="RandomTree.t:S%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:S%s" taxa="@S%s">\n', local_clusters{lc},local_clusters{lc},local_clusters{lc});
                    fprintf(g, '\t\t\t <populationModel id="ConstantPopulation0.t:S%s" spec="ConstantPopulation">\n', local_clusters{lc});
                    fprintf(g, '\t\t\t\t <parameter id="randomPopSize.t:S%s" spec="parameter.RealParameter" name="popSize">0.1</parameter>\n', local_clusters{lc});
                    fprintf(g, '\t\t\t </populationModel>\n');
                    fprintf(g, '\t\t </init>\n');

                end
     %       elseif contains(line, 'insert_mcov')
      %           fprintf(g,'\t\t<covariates id="migration_forward" spec="beast.mascot.glmmodel.Covariate">%s</covariates>\n', sprintf('%f ', mcov_forward));
       %          fprintf(g,'\t\t<covariates id="migration_back" spec="beast.mascot.glmmodel.Covariate">%s</covariates>\n', sprintf('%f ', mcov_back));
            elseif contains(line, 'insert_skyline')
                skyline_shift = [0.1042:2.5/24:2.5];
                
                fprintf(g, strrep(line, 'insert_skyline',num2str(skyline_shift)));
            elseif contains(line, 'insert_rate_shift')
                fprintf(g,'\t\t<rateShifts id="rateShifts.t:S1" spec="mascot.util.InitializedRateShifts" dateFormat="decimal">%s</rateShifts>\n', sprintf('%f ', rate_shift));
                
            elseif contains(line, 'insert_location')
                loc = 'rem';
                for k = 1 :length(dat{6})

                    loc = [loc ',' dat{1}{k} '=' dat{3}{k}];
                end

                loc = strrep(loc, 'rem,','');            
                fprintf(g, strrep(line, 'insert_location',loc));
            elseif contains(line, 'insert_taxa')
                for k = 1 :length(dat{6})
                    fprintf(g, '\t\t\t\t\t\t\t\t<taxon id="%s" spec="beast.evolution.alignment.Taxon"/>\n', dat{1}{k});
                end

     
            elseif contains(line, 'insert_multitrees')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t\t\t\t\t<tree idref="Tree.t:S%s"/>\n', local_clusters{lc});                            
                    fprintf(g, '\t\t\t\t\t\t<rootLength idref="rootLength:S%s"/>\n', local_clusters{lc});
                end
                
% 
%             elseif contains(line, 'insert_priors')
%                 for lc = 1:length(local_clusters)
%                     fprintf(g, '\t\t\t\t<prior id="rootLength.t:S%s" name="distribution" x="@rootLength:S%s">\n', local_clusters{lc}, local_clusters{lc});                                 
%                     fprintf(g, '\t\t\t\t\t<Exponential id="Exponential.S%s" name="distr">\n', local_clusters{lc});        
%                     fprintf(g, '\t\t\t\t\t\t<parameter id="RealParameter.S%s" spec="parameter.RealParameter" lower="0.0" name="mean" upper="0.0">0.01</parameter>\n', local_clusters{lc});
%                     fprintf(g, '\t\t\t\t\t</Exponential>\n');
%                     fprintf(g, '\t\t\t\t</prior>\n');
%                 end                

            elseif contains(line, 'insert_likelihood')
                for lc = 1:length(local_clusters)
                    names = dat{1}(ismember(dat{6},local_clusters{lc}));
                    if length(names)>1
                            fprintf(g, '\t\t\t\t<distribution id="treeLikelihood.S%s" spec="ThreadedTreeLikelihood" data="@S%s" tree="@Tree.t:S%s" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n', local_clusters{lc}, local_clusters{lc}, local_clusters{lc});                

                    end

                end
                



    elseif contains(line, 'insert_operators')
            fprintf(g, '\t\t<operator id="UpDown" spec="UpDownOperator" scaleFactor="0.5" weight="%f">\n', 10);
            for lc = 1:length(local_clusters)
                names = dat{1}(ismember(dat{6},local_clusters{lc}));
                if length(names)>1
                    fprintf(g, '\t\t<up idref="Tree.t:S%s"/>\n', local_clusters{lc});
                end
            end
            fprintf(g, '\t\t</operator>\n');

            for lc = 1:length(local_clusters)
                names = dat{1}(ismember(dat{6},local_clusters{lc}));
                fprintf(g, '\t\t<operator id="rootLengthScaler.t:S%s" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 0.5);
                if length(names)>1
                    fprintf(g, '\t\t<operator id="MascotTreeScaler.t:S%s" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 0.5);
                    fprintf(g, '\t\t<operator id="MascotTreeRootScaler.t:S%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 0.5);
                    if length(names)>2
                        weight = length(names)/length(dat{1})*5;
                        fprintf(g, '\t\t<operator id="MascotUniformOperator.t:S%s" spec="Uniform" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 30*weight);
                        fprintf(g, '\t\t<operator id="MascotSubtreeSlide.t:S%s" spec="SubtreeSlide" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 15*weight);
                        fprintf(g, '\t\t<operator id="MascotNarrow.t:S%s" spec="Exchange" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 15*weight);
                        fprintf(g, '\t\t<operator id="MascotWide.t:S%s" spec="Exchange" isNarrow="false" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 3*weight);
                        fprintf(g, '\t\t<operator id="MascotWilsonBalding.t:S%s" spec="WilsonBalding" tree="@Tree.t:S%s" weight="%f"/>\n', local_clusters{lc}, local_clusters{lc}, 3*weight);
                    end
                end
            end
            
            elseif contains(line, 'insert_logs')
              %  fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');
               % fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');

                 for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t\t<log idref="rootLength:S%s"/>\n', local_clusters{lc});        
                 end

                 for lc = 1 : length(local_clusters)
                    fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:S%s"/>\n',lc,local_clusters{lc});
              %      fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger1:%d" spec="nab.util.MultiTreeStatLogger" heightOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters{lc},local_clusters{lc});
               %     fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger2:%d" spec="nab.util.MultiTreeStatLogger" originOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters{lc},local_clusters{lc});
                 end
             elseif contains(line, 'insert_logtree')

                fprintf(g,'\t\t <logger id="typedTreeloggedsar.t:EBOV_1" spec="Logger" fileName="$(filebase).typed.trees" logEvery="500000" mode="tree">\n');
                fprintf(g,'\t\t\t<log id="structuredTredaslog.t:EBOV_1" spec="nab.multitree.MappedMultitreeMascot" dynamics="@StructuredSkyline.t:S1" immigrationRate="@timeVaryingMigrationRates" multiTreeIntervals="@StructuredMultiTreeIntervals" implementation="indicators">\n');
                fprintf(g,'\t\t\t</log>\n');
                fprintf(g,'\t\t</logger>\n');
                

            else
                fprintf(g, line);
            end
        end
    end
end