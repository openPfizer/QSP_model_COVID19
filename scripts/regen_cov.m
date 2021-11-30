%{ 
    NOTE: Script creates publication ready simulations and figures for the virological endpoints for REGEN-COV nAb cocktail from the REGEN-COV Ph2 and Ph2 trials
    
    Weinreich, D.M.; Sivapalasingam, S.; Norton, T.; Ali, S.; Gao, H.; Bhore, R.; Musser, B.J.; Soo, Y.; Rofail, D.; Im, J., et al. REGN-COV2, a Neutralizing Antibody Cocktail, in Outpatients with Covid-19. N Engl J Med 2021, 384, 238-251, doi:10.1056/NEJMoa2035002.

    Weinreich, D.M.; Sivapalasingam, S.; Norton, T.; Ali, S.; Gao, H.; Bhore, R.; Xiao, J.; Hooper, A.T.; Hamilton, J.D.; Musser, B.J., et al. REGEN-COV Antibody Combination and Outcomes in Outpatients with Covid-19. N Engl J Med 2021, 10.1056/NEJMoa2108163, doi:10.1056/NEJMoa2108163
%}

load("regen_cov.mat")
%% SECTION - Run PBO solutions for REGEN vpop

    data_dictionary_orig = get_data_dictionary(); % load model dictionary
    states_to_store = [1,28];
    statenames_to_store = data_dictionary_orig.species_names(states_to_store,2);
    tmp_parameters = data_dictionary_orig.parameters;
    n_freqav = 12; % sampling frequency of solutions
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arraynomavregn = zeros(n_time_pointsav,length(states_to_store),length(samples_regn)); % initialize solution array
    err_vectornom = [];
    virus_inoculation = 1e1; % initial viral inoculum
    data_dictionary_orig.pknab = 'none';
    parfor ii = 1:length(samples_regn)
        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_regn); % Update data_dictionary for each sample  
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        ind_inf = find(T<=0,1,'last'):find(T<=40,1,'last');
        T = T(ind_inf);
        X = X(ind_inf,:);
        size_sample = length(T(1:n_freqav:end));
        state_arraynomavregn(:,:, ii) = X(1:n_freqav:end,states_to_store);

        % discard if failure and store parameter index
        if size_sample < n_time_pointsav
            err_vectornom = [err_vectornom;ii];
        end
        disp(ii)
    end
    data_dictionary_orig.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
    data_dictionary_orig.parameters.nhill = 1;
    data_dictionary_orig.pknab = 'REGN';
    data_dictionary_orig.nAb0 = 4000;
    [T,X] = function_run_model_noplots(data_dictionary_orig,virus_inoculation,0);
    ind_inf = find(T<=0,1,'last'):find(T<=40*24,1,'last');
    T = T(ind_inf);
    X = X(ind_inf,:);
    T_sample = T(1:n_freqav:end);
    T_sampleav = T_sample;
    % Get Vpeak and T(Vpeak) in hrs
    [maxv,maxiv] = max(squeeze(state_arraynomavregn(:,1,:)));
    tvmax_3s_regn = T_sample(maxiv)/24;
% !SECTION - end Run PBO solutions for REGEN vpop


%% SECTION - Run REGEN-COV nAb simulations for REGEN-COV 8000mg
    udelay = 4; % administration time with respect to exposure
    filter_3s = samples_regn;
    Delayab = repelem(4,length(samples_regn));
    n_vpav = numel(Delayab);
    n_freqav = 12;
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arrayregn = zeros(n_time_pointsav,2,numel(Delayab));
    tmp_parameters = data_dictionary_orig.parameters;
    err_vectorbnab = -1*ones(numel(Delayab),1);
    parfor ii = 1:n_vpav
        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_regn); % Update data_dictionary for each sample
        data_dictionary.parameters.tau_Ab = 24*(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii));
        data_dictionary.pknab = 'REGN';
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (20d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        data_dictionary.nAb0 = 4000;
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        ind_inf = find(T<=0,1,'last'):find(T<=40,1,'last');
        T = T(ind_inf);
        X = X(ind_inf,:);
        T_sample = T(1:n_freqav:end); %% only append 1 to index the size of T instead of all just incase there is anintergration error
        size_sample = length(T(1:n_freqav:end));
        state_arrayregn(:,:, ii) = X(1:n_freqav:end,[1,28]);


        if size_sample < n_time_pointsav
            err_vectorbnab = [err_vectorbnab;ii];
        end
        disp(ii)
    end


    V_nomab = squeeze(state_arraynomavregn(:,1,:));
    V_nomab = repmat(V_nomab,1,numel(Delayab)/length(samples_regn));
    V_treatab = squeeze(state_arrayregn(:,1,:));
    V_nomtab = -1*ones(5,numel(Delayab));
    V_treatabt = -1*ones(5,numel(Delayab));

    IL6_nomab = squeeze(state_arraynomavregn(:,2,:));
    IL6_nomab = repmat(IL6_nomab,1,numel(Delayab)/length(samples_regn));
    IL6_treatab = squeeze(state_arrayregn(:,2,:));
    IL6_nomtab = -1*ones(5,numel(Delayab));
    IL6_treatabt = -1*ones(5,numel(Delayab));
    IL6_nomabmax = -1*ones(1,numel(Delayab));

    for ii = 1:numel(Delayab)
     % Get time vector indices for trial time points
        ind = [find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii)),1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+2,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+4,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+6,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+10,1,'last')];
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
    diff_nombaseregn = log10(V_nomtab)- log10(V_nomtab(1,:));
    vnombase_regn = log10(V_nomtab);
    diff_treatbaseregn = log10(V_treatabt) - log10(V_treatabt(1,:));
    vptreatab_regn = log10(V_treatabt);
    il6nombase_regn = IL6_nomab*mw.il6;
    il6ptreatab_regn = IL6_treatab*mw.il6;
% !SECTION - end Run REGEN-COV nAb simulations for REGEN-COV 8000mg

%% SECTION - REGEN-COV2 viral load time course plots

    figure,
    tl = tiledlayout(3,2);
    nexttile,
    % Baseline Viral load > 1e7
    ind1 = find(vnombase_regn(1,:)>7.0);
    vpnomtab = vnombase_regn(:,ind1);
    vptreatab_regn2 = vptreatab_regn(:,ind1);
    ind2 = 1:length(ind1);

    hold on, 
    xlim([0 8])
    ylim([2 8])
    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')

    % REGEN-COV2 trial data

    regn_pbo_mean = [7.7
                    6.8
                    6.1
                    5.4];
    regn_pbo_err = [7.8
                    7.1
                    6.4
                    5.9];
    regn_8mg_mean = [7.670454545
                    5.795454545
                    4.460227273
                    3.664772727];
    regn_8mg_err = [7.784090909
                    6.193181818
                    4.772727273
                    4.0625];


    errorbar([1,3,5,7],(regn_pbo_mean), regn_pbo_err-regn_pbo_mean,'LineStyle','--','LineWidth',2.5,'Color','k')
    hold on,
    errorbar([0,2,4,6]+1,mean(vpnomtab(1:4,ind2),2),std(vpnomtab(1:4,ind2),[],2)/sqrt(28),'LineWidth',2.5,'Color','k')

    errorbar([1,3,5,7],(regn_8mg_mean), regn_8mg_err-regn_8mg_mean,'LineStyle','--','LineWidth',2.5,'Color',[0,0.4470,0.7410])
    errorbar([0,2,4,6]+1,mean(vptreatab_regn2(1:4,ind2),2),std(vptreatab_regn2(1:4,ind2),[],2)/sqrt(22),'LineWidth',2.5,'Color',[0.8500,0.3250,0.0980])
    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')
    legend('PBO','PBO (model simulation)','Treated 8g REGEN-COV','8g REGEN-COV (model simulation)','Location','SouthWest');
    title('Baseline viral load >10^7')
    ax = gca;
    ax.FontSize = 14;
    grid on

    % Baseline Viral load > 1e6
    ind1 = find(vnombase_regn(1,:)>6.0);
    vpnomtab = vnombase_regn(:,ind1);
    vptreatab_regn2 = vptreatab_regn(:,ind1);
    ind2 = 1:length(ind1);

    nexttile,
    hold on, 
    xlim([0 8])
    ylim([2 8])
    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')

    % REGEN-COV2 trial data

    regn_pbo_mean = [7.4
    6.6
    5.7
    4.95];
    regn_pbo_err = [7.54
    6.88
    6.01
    5.35];
    regn_8mg_mean = [7.4
    5.5
    4.3
    3.3];
    regn_8mg_err = [7.55
    5.95
    4.68
    3.88];


    errorbar([1,3,5,7],(regn_pbo_mean), regn_pbo_err-regn_pbo_mean,'LineStyle','--','LineWidth',2.5,'Color','k')
    errorbar([0,2,4,6]+1,mean(vpnomtab(1:4,ind2),2),std(vpnomtab(1:4,ind2)/sqrt(27),[],2),'LineWidth',2.5,'Color','k')

    hold on,
    errorbar([1,3,5,7],(regn_8mg_mean), regn_8mg_err-regn_8mg_mean,'LineStyle','--','LineWidth',2.5,'Color',[0,0.4470,0.7410])
    errorbar([0,2,4,6]+1,mean(vptreatab_regn2(1:4,ind2),2),std(vptreatab_regn2(1:4,ind2),[],2)/sqrt(34),'LineWidth',2.5,'Color',[0.8500,0.3250,0.0980])

    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')
    title('Baseline viral load >10^6')
    ax = gca;
    ax.FontSize = 14;
    grid on

    % Baseline Viral load > 1e5
    ind1 = find(vnombase_regn(1,:)>5.0);
    vpnomtab = vnombase_regn(:,ind1);
    vptreatab_regn2 = (vptreatab_regn(:,ind1));
    ind2 = 1:length(ind1);

    nexttile,
    hold on, 

    xlim([0 8])
    ylim([2 8])
    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')

    % REGEN-COV2 trial data

    regn_pbo_mean = [6.7
    5.8
    4.7
    4.04];
    regn_pbo_err = [6.89
    6.11
    5.12
    4.4];
    regn_8mg_mean = [6.9
    5.29
    4.01
    3.17];
    regn_8mg_err = [7.07
    5.62
    4.33
    3.53];


    errorbar([1,3,5,7],(regn_pbo_mean), regn_pbo_err-regn_pbo_mean,'LineStyle','--','LineWidth',2.5,'Color','k')
    hold on,
    errorbar([0,2,4,6]+1,mean(vpnomtab(1:4,ind2),2),std(vpnomtab(1:4,ind2)/sqrt(41),[],2),'LineWidth',2,'Color','k')

    errorbar([1,3,5,7],(regn_8mg_mean), regn_8mg_err-regn_8mg_mean,'LineStyle','--','LineWidth',2.5,'Color',[0,0.4470,0.7410])
    errorbar([0,2,4,6]+1,mean(vptreatab_regn2(1:4,ind2),2),std(vptreatab_regn2(1:4,ind2),[],2)/sqrt(45),'LineWidth',2,'Color',[0.8500,0.3250,0.0980])

    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')
    title('Baseline viral load >10^5')
    ax = gca;
    ax.FontSize = 14;
    grid on


    % Baseline viral load > 1e4
    ind1 = find(vnombase_regn(1,:)>4.0);
    vpnomtab = vnombase_regn(:,ind1);
    vptreatab_regn2 = (vptreatab_regn(:,ind1));
    ind2 = 1:length(ind1);


    nexttile,
    hold on, 
    xlim([0 8])
    ylim([2 8])
    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')

    % REGEN-COV2 trial data

    regn_pbo_mean = [6.17
    5.1
    4.3
    3.5];
    regn_pbo_err = [6.4
    5.4
    4.5
    3.8];
    regn_8mg_mean = [6.5
    4.9
    3.8
    2.9];
    regn_8mg_err = [6.7
    5.2
    4.1
    3.26];


    errorbar([1,3,5,7],(regn_pbo_mean), regn_pbo_err-regn_pbo_mean,'LineStyle','--','LineWidth',2.5,'Color','k')
    hold on,
    errorbar([0,2,4,6]+1,mean(vpnomtab(1:4,ind2),2),std(vpnomtab(1:4,ind2)/sqrt(56),[],2),'LineWidth',2.5,'Color','k')

    errorbar([1,3,5,7],(regn_8mg_mean), regn_8mg_err-regn_8mg_mean,'LineStyle','--','LineWidth',2.5,'Color',[0,0.4470,0.7410])
    errorbar([0,2,4,6]+1,mean(vptreatab_regn2(1:4,ind2),2),std(vptreatab_regn2(1:4,ind2),[],2)/sqrt(54),'LineWidth',2.5,'Color',[0.8500,0.3250,0.0980])

    ylabel('Viral Load log_{10} copies/mL')
    xlabel('Days from Start of Treatment')
    title('Baseline viral load >10^4')
    ax = gca;
    ax.FontSize = 14;
    grid on

% !SECTION - end REGEN-COV2 viral load time course plots

%% SECTION - Boostrap Vpop for confidence intervals of REGEN-COV2 nAb trial
    maxil6ptreat_regn = max(il6ptreatab_regn);
    maxil6nombase_regn = max(il6nombase_regn);
    n_bs = 1000;
    rrr_regn = -1*ones(1000,1);
    pbo_regn = -1*ones(1000,1);
    trt_regn = -1*ones(1000,1);

    for ii = 1:n_bs
        bsind = randsample(1:length(samples_regn),length(samples_regn),true);
        event_pbo = find(maxil6nombase_regn(bsind)>=40);
        event_trt = find(maxil6ptreat_regn(bsind)>=40);
        rrr_regn(ii) = 1 - length(event_trt)/length(event_pbo);
        pbo_regn(ii) = length(event_pbo)/length(samples_regn);
        trt_regn(ii) = length(event_trt)/length(samples_regn);
    end
% !SECTION - Boostrap Vpop for confidence intervals of nAb trial
    mean(rrr_regn)
    max(rrr_regn)
    min(rrr_regn)

%% SECTION - Calculate and plot delta change from baseline for subgroups from REGEN-COV2 trial

        % Total population - 1e4 
            ind1 = find(vnombase_regn(1,:)>4.0);
            vpnomtab = vnombase_regn(:,ind1);
            vptreatab_regn2 = (vptreatab_regn(:,ind1));
            ind2 = 1:length(ind1);

            pbo_vregn = vpnomtab(3,:);
            pbo_vregn(pbo_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            trt_vregn = vptreatab_regn2(3,:);
            trt_vregn(trt_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVregn_bs = -1*ones(n_bs,size(pbo_vregn,1));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(ind2),length(ind2),true);
                pbo_bs = pbo_vregn(:,ind_bs);
                trt_bs = trt_vregn(:,ind_bs);
                diff_meanVregn_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
            end

            diff_meanVregn.x1e4 = diff_meanVregn_bs;

        % Total population - 1e5 
            ind1 = find(vnombase_regn(1,:)>5.0);
            vpnomtab = vnombase_regn(:,ind1);
            vptreatab_regn2 = (vptreatab_regn(:,ind1));
            ind2 = 1:length(ind1);

            pbo_vregn = vpnomtab(3,:);
            pbo_vregn(pbo_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            trt_vregn = vptreatab_regn2(3,:);
            trt_vregn(trt_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVregn_bs = -1*ones(n_bs,size(pbo_vregn,1));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(ind2),length(ind2),true);
                pbo_bs = pbo_vregn(:,ind_bs);
                trt_bs = trt_vregn(:,ind_bs);
                diff_meanVregn_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
            end

            diff_meanVregn.x1e5 = diff_meanVregn_bs;

        % Total population - 1e6
            ind1 = find(vnombase_regn(1,:)>6.0);
            vpnomtab = vnombase_regn(:,ind1);
            vptreatab_regn2 = (vptreatab_regn(:,ind1));
            ind2 = 1:length(ind1);
            
            pbo_vregn = vpnomtab(3,:);
            pbo_vregn(pbo_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            trt_vregn = vptreatab_regn2(3,:);
            trt_vregn(trt_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVregn_bs = -1*ones(n_bs,size(pbo_vregn,1));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(ind2),length(ind2),true);
                pbo_bs = pbo_vregn(:,ind_bs);
                trt_bs = trt_vregn(:,ind_bs);
                diff_meanVregn_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
            end

            diff_meanVregn.x1e6 = diff_meanVregn_bs;

        % Total population - 1e7
            ind1 = find(vnombase_regn(1,:)>7.0);
            vpnomtab = vnombase_regn(:,ind1);
            vptreatab_regn2 = (vptreatab_regn(:,ind1));
            ind2 = 1:length(ind1);
            
            pbo_vregn = vpnomtab(3,:);
            pbo_vregn(pbo_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            trt_vregn = vptreatab_regn2(3,:);
            trt_vregn(trt_vregn<=1e0) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVregn_bs = -1*ones(n_bs,size(pbo_vregn,1));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(ind2),length(ind2),true);
                pbo_bs = pbo_vregn(:,ind_bs);
                trt_bs = trt_vregn(:,ind_bs);
                diff_meanVregn_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
            end

            diff_meanVregn.x1e7 = diff_meanVregn_bs;


        % >1e4
        regn_pbo_mean = [6.17
                        5.1
                        4.3
                        3.5];
        regn_8mg_mean = [6.5
                        4.9
                        3.8
                        2.9];
        diff_regn_1e4 = (regn_8mg_mean(1)-regn_8mg_mean(3)) - (regn_pbo_mean(1) - regn_pbo_mean(3));


        % > 1e5
        regn_pbo_mean = [6.7
                        5.8
                        4.7
                        4.04];
        regn_8mg_mean = [6.9
                        5.29
                        4.01
                        3.17];
        diff_regn_1e5 = (regn_8mg_mean(1)-regn_8mg_mean(3)) - (regn_pbo_mean(1) - regn_pbo_mean(3));

        % > 1e6
        regn_pbo_mean = [7.4
                        6.6
                        5.7
                        4.95];
        
        regn_8mg_mean = [7.4
                        5.5
                        4.3
                        3.3];
        diff_regn_1e6 = (regn_8mg_mean(1)-regn_8mg_mean(3)) - (regn_pbo_mean(1) - regn_pbo_mean(3));

        % > 1e7
        regn_pbo_mean = [7.7
                    6.8
                    6.1
                    5.4];
        regn_8mg_mean = [7.67
                    5.8
                    4.46
                    3.66];
        diff_regn_1e7 = (regn_8mg_mean(1)-regn_8mg_mean(3)) - (regn_pbo_mean(1) - regn_pbo_mean(3));


        nexttile,
        plt(1) = scatter(mean(diff_meanVregn.x1e4),diff_regn_1e4,512,[0.8500,0.3250,0.0980],'filled','Marker','d');
        hold on,
        plt(2) = scatter(mean(diff_meanVregn.x1e5),diff_regn_1e5,512,[0.8500,0.3250,0.0980],'filled','Marker','o');
        hold on,
        plt(3) = scatter(mean(diff_meanVregn.x1e6),diff_regn_1e6,512,[0.8500,0.3250,0.0980],'filled','Marker','^');
        hold on,
        plt(4) = scatter(mean(diff_meanVregn.x1e7),diff_regn_1e7,512,[0.8500,0.3250,0.0980],'filled','Marker','v');
        hold on,
        % scatter([mean(diff_meanVregn.x1e4),mean(diff_meanVregn.x1e5),mean(diff_meanVregn.x1e6),mean(diff_meanVregn.x1e7)],[diff_regn_1e4,diff_regn_1e5,diff_regn_1e6,diff_regn_1e7],512,0*[0,0.4470,0.7410; 0.8500,0.3250,0.0980; 0.9290,0.6940,0.1250; 0.4940,0.1840,0.5560],'filled','Marker','d')
        e2 = errorbar(mean(diff_meanVregn.x1e4),diff_regn_1e4, [],[],(mean(diff_meanVregn.x1e4) - prctile(diff_meanVregn.x1e4,0.5)),(prctile(diff_meanVregn.x1e4,99.5) - mean(diff_meanVregn.x1e4)),'LineWidth',2,'LineStyle','none','Color','k');
        errorbar(mean(diff_meanVregn.x1e5),diff_regn_1e5, [],[],(mean(diff_meanVregn.x1e5) - prctile(diff_meanVregn.x1e5,0.5)),(prctile(diff_meanVregn.x1e5,99.5) - mean(diff_meanVregn.x1e5)),'LineWidth',2,'LineStyle','none','Color','k') 
        errorbar(mean(diff_meanVregn.x1e6),diff_regn_1e6, [],[],(mean(diff_meanVregn.x1e6) - prctile(diff_meanVregn.x1e6,0.5)),(prctile(diff_meanVregn.x1e6,99.5) - mean(diff_meanVregn.x1e6)),'LineWidth',2,'LineStyle','none','Color','k')
        errorbar(mean(diff_meanVregn.x1e7),diff_regn_1e7, [],[],(mean(diff_meanVregn.x1e7) - prctile(diff_meanVregn.x1e7,0.5)),(prctile(diff_meanVregn.x1e7,99.5) - mean(diff_meanVregn.x1e7)),'LineWidth',2,'LineStyle','none','Color','k') 
        xticks(0:0.25:2)
        yticks(0:0.25:2)
        xlim([0 2.0])
        ylim([0 2.0])
        ylabel({'REGEN-COV Ph2 \Delta V,',' -log_{10} RNA copies/mL at Day 5'})
        xlabel('Simulated \Delta V, -log_{10} RNA copies/mL at Day 5')
        s1 = scatter(-10,-10,100,0*[0,0.4470,0.7410],'filled','Marker','d');
        s2 = scatter(-10,-10,100,0*[0.8500,0.3250,0.0980],'filled','Marker','o');
        s3 = scatter(-10,-10,100,0*[0.9290,0.6940,0.1250],'filled','Marker','s');
        s4 = scatter(-10,-10,100,0*[0.4940,0.1840,0.5560],'filled','Marker','+');
        grid on
        set(gca,'FontSize',14)
        rf = refline(1,0);
        rf.LineWidth = 2;
        rf.LineStyle = '--';
        rf.Color = 'k';
        [leg] = legend([plt(1),plt(2),plt(3),plt(4)],'>4log_{10}','>5log_{10}','>6log_{10}','>7log_{10}','Location','SouthEast','FontSize',12);
        leg.Title.String = {"Baseline viral load", "(RNA copies/mL)"};
        legend boxoff

% !SECTION - Calculate and plot delta change from baseline for subgroups


%% SECTION - Run REGEN-COV nAb simulations for REGEN-COV 2400mg
    udelay = 4; % administration time with respect to exposure
    filter_3s = samples_regn;
    Delayab = repelem(4,length(samples_regn));
    n_vpav = numel(Delayab);
    n_freqav = 12;
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arrayregn3 = zeros(n_time_pointsav,2,numel(Delayab));
    tmp_parameters = data_dictionary_orig.parameters;
    err_vectorbnab = -1*ones(numel(Delayab),1);
    parfor ii = 1:n_vpav
        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_regn); % Update data_dictionary for each sample
        data_dictionary.parameters.tau_Ab = 24*(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii));
        data_dictionary.pknab = 'REGN';
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (20d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        data_dictionary.nAb0 = 1200;
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        ind_inf = find(T<=0,1,'last'):find(T<=40,1,'last');
        T = T(ind_inf);
        X = X(ind_inf,:);
        T_sample = T(1:n_freqav:end); %% only append 1 to index the size of T instead of all just incase there is anintergration error
        size_sample = length(T(1:n_freqav:end));
        state_arrayregn3(:,:, ii) = X(1:n_freqav:end,[1,28]);


        if size_sample < n_time_pointsav
            err_vectorbnab = [err_vectorbnab;ii];
        end
        disp(ii)
    end


    V_nomab = squeeze(state_arraynomavregn(:,1,:));
    V_nomab = repmat(V_nomab,1,numel(Delayab)/length(samples_regn));
    V_treatab = squeeze(state_arrayregn3(:,1,:));
    V_nomtab = -1*ones(5,numel(Delayab));
    V_treatabt = -1*ones(5,numel(Delayab));

    IL6_nomab = squeeze(state_arraynomavregn(:,2,:));
    IL6_nomab = repmat(IL6_nomab,1,numel(Delayab)/length(samples_regn));
    IL6_treatab = squeeze(state_arrayregn3(:,2,:));
    IL6_nomtab = -1*ones(5,numel(Delayab));
    IL6_treatabt = -1*ones(5,numel(Delayab));
    IL6_nomabmax = -1*ones(1,numel(Delayab));

    for ii = 1:numel(Delayab)
     % Get time vector indices for trial time points
        ind = [find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii)),1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+2,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+4,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+6,1,'last');
        find(T_sampleav/24<=(tvmax_3s_regn(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delayab(ii))+10,1,'last')];
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
    diff_nombaseregn = log10(V_nomtab)- log10(V_nomtab(1,:));
    vnombase_regn3 = log10(V_nomtab);
    diff_treatbaseregn = log10(V_treatabt) - log10(V_treatabt(1,:));
    vptreatab_regn3 = log10(V_treatabt);
    il6nombase_regn = IL6_nomab*mw.il6;
    il6ptreatab_regn = IL6_treatab*mw.il6;
% !SECTION - end Run REGEN-COV nAb simulations for REGEN-COV 2400mg

%% SECTION - Boostrap Vpop for confidence intervals of REGEN-COV Ph3 nAb trial
    maxil6ptreat_regn = max(il6ptreatab_regn);
    maxil6nombase_regn = max(il6nombase_regn);
    n_bs = 1000;
    rrr_regn = -1*ones(1000,1);
    pbo_regn = -1*ones(1000,1);
    trt_regn = -1*ones(1000,1);

    for ii = 1:n_bs
        bsind = randsample(1:length(samples_regn),length(samples_regn),true);
        event_pbo = find(maxil6nombase_regn(bsind)>=40);
        event_trt = find(maxil6ptreat_regn(bsind)>=40);
        rrr_regn(ii) = 1 - length(event_trt)/length(event_pbo);
        pbo_regn(ii) = length(event_pbo)/length(samples_regn);
        trt_regn(ii) = length(event_trt)/length(samples_regn);
    end

    mean(rrr_regn)
    max(rrr_regn)
    min(rrr_regn)

    % Total population - 1e4 

    pbo_vregn = vnombase_regn(4,:);
    trt_vregn = vptreatab_regn3(4,:);
    n_bs = 1000;
    diff_meanVregn_bs = -1*ones(n_bs,size(pbo_vregn,1));
    for kk = 1:n_bs
        ind_bs = randsample(1:length(samples_regn),length(samples_regn),true);
        pbo_bs = pbo_vregn(:,ind_bs);
        trt_bs = trt_vregn(:,ind_bs);
        diff_meanVregn_bs(kk,:) = (mean((pbo_bs),2) - mean((trt_bs),2))';
    end

    diff_meanVregn3.x1e4 = diff_meanVregn_bs;
    mean(diff_meanVregn3.x1e4)
% !SECTION - Boostrap Vpop for confidence intervals of REGEN-COV Ph3 nAb trial
    nexttile,
    scatter(1,0.86,512,'filled','Marker','d')
    hold on
    vsplot = violinplot([-1;diff_meanVregn3.x1e4],[0;ones(size(diff_meanVregn3.x1e4))],'ViolinAlpha',1,'ViolinColor',[0.8500,0.3250,0.0980],'ShowData',false);
    e2 = errorbar(1,0.86, 0.72-0.86,1-0.86,'LineWidth',1.5,'LineStyle','none','Color','k');  
    grid on
    xlabel('[Treatment starts at Day 1]')
    ylabel('log_{10}[reduction in viral load] at Day 7')
    ylim([0.5 1.25])
    xlim([0.75 2.5])
    set(gca,'FontSize',16)
    grid on
    row1 = {'2.4g REGEN-COV','Simulated'};
    row2 = {'Ph 3 trial', '2.4g REGEN-COV'};
    labelArray = [row1; row2]; 
    labelArray = strjust(pad(labelArray),'center');
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    ax = gca(); 
    ax.XTickLabel = tickLabels;
    set(gcf, 'Position',  [300   527   1100   1000])
