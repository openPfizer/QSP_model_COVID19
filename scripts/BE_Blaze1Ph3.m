%{ 
    NOTE: Script creates publication ready simulations and figures for the virological endpoints for 2800mg bamlanivimab + 2800mg etesevimab nAb cocktail from the Blaze-1 Ph3 clinical trial 
    
    [Dougan, M.; Nirula, A.; Azizad, M.; Mocherla, B.; Gottlieb, R.L.; Chen, P.; Hebert, C.; Perry, R.; Boscia, J.; Heller, B., et al. Bamlanivimab plus Etesevimab in Mild or Moderate Covid-19. N Engl J Med 2021, 385, 1382-1392, doi:10.1056/NEJMoa2102685]

    Emergency Use Authorization (EUA) for Bamlanivimab 700 mg and Etesevimab 1400 mg IV Administered Together Center for Drug Evaluation and Research (CDER) Review. 131.

%}
load("blaze1.mat")
%% SECTION Run PBO solutions for Blaze-1 Ph3 vpop
    data_dictionary_orig = get_data_dictionary(); % load model dictionary
    states_to_store = [1,28];
    statenames_to_store = data_dictionary_orig.species_names(states_to_store,2);
    tmp_parameters = data_dictionary_orig.parameters;
    n_freqav = 12; % sampling frequency of solutions
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arraynomav = zeros(n_time_pointsav,length(states_to_store),length(samples_indmin)); % initialize solution array
    err_vectornom = [];
    virus_inoculation = 1e1; % initial viral inoculum
    data_dictionary_orig.pknab = 'none';
    parfor ii = 1:length(samples_indmin)
        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_indmin); % Update data_dictionary for each sample  
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        ind_inf = find(T<=0,1,'last'):find(T<=40,1,'last');
        T = T(ind_inf);
        X = X(ind_inf,:);
        size_sample = length(T(1:n_freqav:end));
        state_arraynomav(:,:, ii) = X(1:n_freqav:end,states_to_store);

        % discard if failure and store parameter index
        if size_sample < n_time_pointsav
            err_vectornom = [err_vectornom;ii];
        end
        disp(ii)
    end
    data_dictionary_orig.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
    data_dictionary_orig.parameters.nhill = 1;
    data_dictionary_orig.pknab = 'LY';
    data_dictionary_orig.nAb0 = 2800;
    [T,X] = function_run_model_noplots(data_dictionary_orig,virus_inoculation,0);
    ind_inf = find(T<=0,1,'last'):find(T<=40*24,1,'last');
    T = T(ind_inf);
    X = X(ind_inf,:);
    T_sample = T(1:n_freqav:end);
    T_sampleav = T_sample;
    % Get Vpeak and T(Vpeak) in hrs
    [maxv,maxiv] = max(squeeze(state_arraynomav(:,1,:)));
    tvmax_3s = T_sample(maxiv)/24;
% !SECTION - end Run PBO solutions for Blaze-1 Ph3 vpop


%% SECTION - Run Blaze-1 Ph3 nAb solutions
    % NOTE - Only simulating nAb administration 4d post viral load peak; comparable to Blaze-1 Ph3 time from symptom onset
    data_dictionary_orig = get_data_dictionary(); % load model dictionary
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arraybnab = zeros(n_time_pointsav,2,length(samples_indmin));
    err_vectorbnab = [];
    Delayab = repelem(4,length(samples_indmin));
    tmp_parameters = data_dictionary_orig.parameters;
    filter_3s = samples_indmin;
    parfor ii = 1:length(samples_indmin)

        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_indmin); % Update data_dictionary for each sample
        data_dictionary.parameters.tau_Ab = 24*(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii));
        data_dictionary.pknab = 'LY';
        data_dictionary.nAb0 = 2800;
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (20d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        ind_inf = find(T<=0,1,'last'):find(T<=40,1,'last');
        T = T(ind_inf);
        X = X(ind_inf,:);
        T_sample = T(1:n_freqav:end); %% only append 1 to index the size of T instead of all just incase there is anintergration error
        size_sample = length(T(1:n_freqav:end));
        state_arraybnab(:,:, ii) = X(1:n_freqav:end,[1,28]);


        if size_sample < n_time_pointsav
            err_vectorbnab = [err_vectorbnab;ii];
        end
        disp(ii)
    end


    mw = data_dictionary_orig.mw;

    V_nomab = squeeze(state_arraynomav(:,1,:));
    V_nomab = repmat(V_nomab,1,numel(Delayab)/length(samples_indmin));
    V_treatab = squeeze(state_arraybnab(:,1,:));
    V_nomtab = -1*ones(5,numel(Delayab));
    V_treatabt = -1*ones(5,numel(Delayab));

    IL6_nomab = squeeze(state_arraynomav(:,2,:));
    IL6_nomab = repmat(IL6_nomab,1,numel(Delayab)/length(samples_indmin));
    IL6_treatab = squeeze(state_arraybnab(:,2,:));
    IL6_nomtab = -1*ones(5,numel(Delayab));
    IL6_treatabt = -1*ones(5,numel(Delayab));
    IL6_nomabmax = -1*ones(1,numel(Delayab));

    for ii = 1:numel(Delayab)
     % Get time vector indices for trial time points
        ind = [find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii)),1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+2,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+4,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+6,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+10,1,'last')];
        V_nomb = V_nomab(ind,ii);
        V_nomb(V_nomb<1e0) = 1;
        V_nomtab(:,ii) = V_nomb;

        V_treatabb = V_treatab(ind,ii);
        V_treatabb(V_treatabb<1e0) = 1;
        V_treatabt(:,ii) = V_treatabb;

        IL6_nomb = IL6_nomab(ind,ii);
        IL6_nomtab(:,ii) = IL6_nomb;

        IL6_treatabb = IL6_treatab(ind,ii);
        IL6_treatabt(:,ii) = IL6_treatabb;
    end
    diff_nombase = log10(V_nomtab)- log10(V_nomtab(1,:));
    vnombase = log10(V_nomtab);
    diff_treatbase = log10(V_treatabt) - log10(V_treatabt(1,:));
    vptreatab = log10(V_treatabt);
    il6nombase = IL6_nomab*mw.il6;
    il6ptreatab = IL6_treatab*mw.il6;
% !SECTION - end Run Blaze-1 Ph3 nAb solutions

%% SECTION - Boostrap Vpop for confidence intervals of Blaze-1 Ph3 trial
    maxil6ptreatab = max(il6ptreatab);
    maxil6nombase = max(il6nombase);
    n_bs = 1000;
    rrr_ly = -1*ones(1000,1);
    pbo_ly = -1*ones(1000,1);
    trt_ly = -1*ones(1000,1);

    for ii = 1:n_bs
        bsind = randsample(1:length(samples_indmin),length(samples_indmin),true);
        event_pbo = find(maxil6nombase(bsind)>=40);
        event_trt = find(maxil6ptreatab(bsind)>=40);
        rrr_ly(ii) = 1 - length(event_trt)/length(event_pbo);
        pbo_ly(ii) = length(event_pbo)/length(samples_indmin);
        trt_ly(ii) = length(event_trt)/length(samples_indmin);
    end

    pbo_vly = round(vnombase([3,4],:),2);
    trt_vly = round(vptreatab([3,4],:),2);
    n_bs = 5000;
    diff_meanVly_bs = -1*ones(n_bs,size(pbo_vly,1));
    for kk = 1:n_bs
        ind_bs = randsample(1:length(samples_indmin),length(samples_indmin),true);
        pbo_bs = pbo_vly(:,ind_bs);
        trt_bs = trt_vly(:,ind_bs);
        diff_meanVly_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
    end

    median(rrr_ly)
    max(rrr_ly)
    min(rrr_ly)
    mean(diff_meanVly_bs)
% !SECTION - end Boostrap Vpop for confidence intervals of Blaze-1 Ph3 trial

%% SECTION - Plot mean viral load time course and viral load reduction at Day 7
    figure, 
    tl = tiledlayout(1,2);
        nexttile,
        expt = [1, 3, 5, 7, 11];
        pl1 = plot(expt, [6.52, 5.74, 4.68, 4.05, 2.69],'LineStyle','--','Marker','+','Color','k','LineWidth',1.5,'MarkerSize',16);
        hold on,
        plot([0,2,4,6,10]+1,mean(vnombase,2),'LineWidth',2,'Color','k')
        plot(expt, [6.51, 5.04, 3.85, 2.87, 2.21],'LineStyle','--','Marker','+','Color',[0,0.4470,0.7410],'LineWidth',1.5,'MarkerSize',16)
        plot([0,2,4,6,10]+1,mean(vptreatab,2),'LineWidth',2,'Color',[0.8500,0.3250,0.0980])
        ylim([1 8])
        xlim([0 12])
        xticks([1 3 5 7 11])
        ylabel('Viral load log_{10} RNA copies/mL')
        xlabel('Days from Start of Treatment')
        hold on
        ax = gca;
        ax.FontSize = 14;
        grid on
        lg = legend('PBO Blaze-1 Ph3','PBO (mean virtual population)','Treated Blaze-1 Ph3','Treated (mean virtual population)');
        lg.Title.String = '4d post symptom onset';

    nexttile,
        scatter(1,1.20,512,'filled','Marker','d')
        hold on,
        vsplot = violinplot([-1*ones(size(diff_meanVly_bs(:,2)));diff_meanVly_bs(:,2)], [zeros(size(diff_meanVly_bs(:,2)));ones(size(diff_meanVly_bs(:,2)))],'ViolinAlpha',1,'ViolinColor',[0.8500,0.3250,0.0980],'BoxColor',[0,0,0],'ShowData',false);
        e1 = errorbar(1,1.20, (1.20 - 0.94), (1.46 - 1.20),'LineWidth',3,'LineStyle','none','Color','k');   
        grid on
        xlabel('[Treatment starts at Day 1]')
        ylabel('log_{10}[reduction in viral load] at Day 7')
        ylim([0.8 1.75])
        xlim([0.5 2.5])
        set(gca,'FontSize',16)
        grid on
        row1 = {'Blaze-1','Model'};
        row2 = {'Ph3', 'simulation'};
        labelArray = [row1; row2]; 
        labelArray = strjust(pad(labelArray),'center');
        tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
        ax = gca(); 
        ax.XTickLabel = tickLabels;
        set(gcf, 'Position',  [300   527   1200   500])

% !SECTION - end - Plot mean viral load time course and viral load reduction at Day 7