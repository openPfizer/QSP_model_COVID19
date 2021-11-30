%{ 
    NOTE: Script creates publication ready simulations and figures for the virological endpoints for 800mg molnupiravir 5d BID treatment from the Ph2a trial 

    Virological endpoints extracted from Fischer, W.; Eron, J.J.; Holman, W.; Cohen, M.S.; Fang, L.; Szewczyk, L.J.; Sheahan, T.P.; Baric, R.; Mollan, K.R.; Wolfe, C.R., et al. Molnupiravir, an Oral Antiviral Treatment for COVID-19. medRxiv 2021, 10.1101/2021.06.17.21258639, doi:10.1101/2021.06.17.21258639

    Pharmacokinetics for EIDD-1931 extracted from Painter, W.P.; Holman, W.; Bush, J.A.; Almazedi, F.; Malik, H.; Eraut, N.; Morin, M.J.; Szewczyk, L.J.; Painter, G.R. Human Safety, Tolerability, and Pharmacokinetics of Molnupiravir, a Novel Broad-Spectrum Oral Antiviral Agent with Activity Against SARS-CoV-2. Antimicrob Agents Chemother 2021, 10.1128/AAC.02428-20, doi:10.1128/AAC.02428-20.
    
    %}
    
load("eidd.mat")

%% SECTION - Load and format molnupiravir PK

    filename = "molnupiravir_PK.csv";
    molnupiravir_pk = readtable(filename,'TextType','string','TreatAsEmpty',"NA");

    filename = "molnupiravir_baselineVL.csv";
    molnupiravir_bsl = readtable(filename,'TextType','string','TreatAsEmpty',"NA");

    filename = "molnupiravir_VLfrombaseline.csv";
    molnupiravir_dbsl = readtable(filename,'TextType','string','TreatAsEmpty',"NA");

    molnupiravir_vl = table('Size',[5,5],'VariableTypes',repmat("double",1,5),'VariableNames',molnupiravir_dbsl.Properties.VariableNames(1:end));
    molnupiravir_vl.Time = molnupiravir_dbsl.Time;
    for ii = 2:length(molnupiravir_dbsl.Properties.VariableNames)
        molnupiravir_vl.(molnupiravir_dbsl.Properties.VariableNames{ii}) = molnupiravir_bsl.(molnupiravir_dbsl.Properties.VariableNames{ii}) + molnupiravir_dbsl.(molnupiravir_dbsl.Properties.VariableNames{ii});
    end

% !SECTION - end Load and format molnupiravir PK

%% SECTION - Run PBO solutions for molnupiravir virtual population
    
    data_dictionary_orig = get_data_dictionary(); % load model dictionary
    tmp_parameters = data_dictionary_orig.parameters;
    n_freqav = 12; % sampling frequency of solutions
    n_time_pointsav = length(0:0.1:(100*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arraynom_merck = zeros(n_time_pointsav,2,length(samples_merck)); % initialize solution array
    err_vectornom = [];
    virus_inoculation = 1e1; % initial viral inoculum
    parfor ii = 1:length(samples_merck)
         [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_merck); % Update data_dictionary for each sample  
        data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
        data_dictionary.parameters.nhill = 1;
        data_dictionary.pknab = 'none';
        [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
        T = T/24;
        size_sample = length(T(1:n_freqav:end));
        state_arraynom_merck(:,:, ii) = X(1:n_freqav:end,[1,28]);

        % discard if failure and store parameter index
        if size_sample < n_time_pointsav
            err_vectornom = [err_vectornom;ii];
        end
        disp(ii)
    end
    data_dictionary_orig.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d = 24*21hr)
    data_dictionary_orig.parameters.nhill = 1;
    data_dictionary_orig.pknab = 'none';
    [T,X] = function_run_model_noplots(data_dictionary_orig,virus_inoculation,0);
    T_sample = T(1:n_freqav:end);
    T_sampleav = T_sample;
    % Get Vpeak and T(Vpeak) in hrs
    [maxv,maxiv] = max(squeeze(state_arraynom_merck(:,1,:)));
    tvmax_3s_merck = T_sample(maxiv)/24;
% !SECTION - end Run PBO solutions for molnupiravir virtual population

%% SECTION - Run molnupiravir simulations

    % SECTION - Load clinical PK approximation

        cmap = [  0.8431    0.0980    0.1098
        0.9290    0.6940    0.1250
        0.1686    0.5137    0.7294];
        figure,
        plot(molnupiravir_pk.Time, molnupiravir_pk.x800mg,'LineWidth',2)
        hold on, 
        plot(molnupiravir_pk.Time, molnupiravir_pk.x600mg,'LineWidth',2)
        plot(molnupiravir_pk.Time, molnupiravir_pk.x400mg,'LineWidth',2)
        plot(molnupiravir_pk.Time, molnupiravir_pk.x300mg,'LineWidth',2)
        plot(molnupiravir_pk.Time, molnupiravir_pk.x200mg,'LineWidth',2)
        plot(molnupiravir_pk.Time, molnupiravir_pk.x100mg,'LineWidth',2)
        plot(molnupiravir_pk.Time, molnupiravir_pk.x50mg,'LineWidth',2)
        ylabel('EIDD-1931 plasma concentrations [ng/mL]')
        xlabel('Time (hours)')
        grid on
        set(gca,'FontSize',14)
        ylim([0 1.25*max(molnupiravir_pk.x800mg)])
        legend('800mg','400mg','300mg','200mg','100mg','50mg')

        % Create molnupiravir PK array to pass to ODE
        molnupiravir_pk = [molnupiravir_pk(1,:);repmat(molnupiravir_pk(2:end,:),10,1)];
        molnupiravir_pk.Time = ([0,repelem(0:9,length(unique(molnupiravir_pk.Time))-1).*(12*ones(1,10*length(unique(molnupiravir_pk.Time(2:end)))))])' + molnupiravir_pk.Time;

    % !SECTION - end Load clinical PK approximation

    fields_to_run = {'x800mg'}; % dose regimens to run
    filter_3s = samples_merck;
    % NOTE - Current assumption Imax = 1
    data_dictionary_orig = get_data_dictionary(); % load model dictionary
    tmp_parameters = data_dictionary_orig.parameters;
    virus_inoculation = 1e1; % initial viral inoculum
    for kk = fields_to_run'
        udelay = 4; % administration time with respect to viral peak
        uimax = 0.999; % Imax
        delay = repmat(repelem(udelay,1,size(filter_3s,1)),1,numel(uimax)); % array of delays for each sample
        imax = repmat(repelem(uimax,1,size(filter_3s,1)),numel(udelay),1); % array of imax for each sample
        Delay = delay';
        Imax = imax';
        n_vpav = numel(Delay);
        n_freqav = 12;
        n_time_pointsav = length(0:0.1:(100*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency- every 50 points
        state_arrayav_mk = NaN(n_time_pointsav,2,numel(Delay)); % initialize solution array for anti-viral simulations
        err_vectorav = -1*ones(numel(Delay),1);
        parfor ii = 1:n_vpav
            [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,ii,perturbation_name_vector,samples_merck);
            data_dictionary.parameters.Imax = Imax(ii);
            data_dictionary.parameters.delay = (tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii));
            data_dictionary.parameters.avduration = 12;
            data_dictionary.parameters.IC50_av = 1086; % NOTE - Estimated molnupiravir IC50
            data_dictionary.parameters.nhill = 1;
            data_dictionary.pknab = 'none';
            data_dictionary.pkcon = molnupiravir_pk.(kk{1}); 
            data_dictionary.pkcont = molnupiravir_pk.Time; %convert PK time from days to hr
            data_dictionary.pkav = true; % settings for using mean clinical EIDD-1931 PK
            data_dictionary.parameters.tau_Ab_Igx = 24*21; % setting for when endogenous Ab response turns on (21d post infection = 24*21hr) 

            virus_inoculation = 1e1; % initial viral inoculum
            [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,0);
            T = T/24;
            size_sample = length(T(1:n_freqav:end));
            
            disp(ii)

            if size_sample < n_time_pointsav
                err_vectorav(ii) = ii; 
            else
                state_arrayav_mk(:,:, ii) = X(1:n_freqav:end,[1,28]);
            end
        end
        simav_merck.(kk{1}) = state_arrayav_mk;
        disp(kk{1})
    end
% !SECTION - end Run molnupiravir simulations

%% SECTION - Analyze solutions/Compare PBO and TRT
    for mm = fieldnames(simav_merck)'
        % Initialize viral Load arrays
        mw = data_dictionary_orig.mw; % molecular weights for species
        V_nom_merck = squeeze(state_arraynom_merck(:,1,:)); % get viral load
        V_nom_merck = repmat(V_nom_merck,1,numel(Delay)/size(filter_3s,1)); % adjust to same dimensions as anti-viral solutions
        V_treatav_merck = squeeze(simav_merck.(mm{1})(:,1,:)); % treated viral load
        V_nomtav_merck = -1*ones(7,numel(Delay)); % initialize array to store solutions at trial time points
        V_treatavt_merck = -1*ones(7,numel(Delay));

        % Initialize IL-6 arrays
        IL6_nomav_merck = squeeze(state_arraynom_merck(:,2,:));
        IL6_nomav_merck = repmat(IL6_nomav_merck,1,numel(Delay)/size(filter_3s,1));
        IL6_treatav_merck = squeeze(simav_merck.(mm{1})(:,2,:));
        IL6_nomtav_merck = -1*ones(7,numel(Delay));
        IL6_treatavt_merck = -1*ones(7,numel(Delay));
        IL6_nomavmax_merck = -1*ones(1,numel(Delay));

        % Get viral load and IL-6 at trial time points
        for ii = 1:numel(Delay)
            % Get time vector indices for trial time points
            ind = [find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii)),1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+2,1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+3,1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+4,1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+5,1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+6,1,'last');
            find(T_sampleav/24<=(tvmax_3s_merck(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+10,1,'last')];
            V_nomb = V_nom_merck(ind,ii);
            V_nomb(V_nomb<=1e0) = 1; % in order to set log(V)<0 = 0
            V_nomtav_merck(:,ii) = V_nomb;

            V_treatavb = V_treatav_merck(ind,ii);
            V_treatavb(V_treatavb<=1e0) = 1; % in order to set log(V)<0 = 0
            V_treatavt_merck(:,ii) = V_treatavb;

            IL6_nomb = IL6_nomav_merck(ind,ii);
            IL6_nomtav_merck(:,ii) = IL6_nomb;

            IL6_treatavb = IL6_treatav_merck(ind,ii);
            IL6_treatavt_merck(:,ii) = IL6_treatavb;
            
        end

        Vtime.(mm{1}).V_nomtav_merck = V_nomtav_merck;
        Vtime.(mm{1}).V_treatavt_merck = V_treatavt_merck;
        Vtime.(mm{1}).V_nom_merck = V_nom_merck;
        Vtime.(mm{1}).V_treatav_merck = V_treatav_merck;

        % Calculate differences in viral load and IL-6 for PBO_i vs TRT_i (paired)
        udelay = unique(delay);
        uimax = unique(imax);
        diff_Vav = cell(length(uimax),length(udelay));
        diff_IL6av = cell(length(uimax),length(udelay));
        diff_IL6maxav = cell(length(uimax),length(udelay));
        IL6maxav = cell(length(uimax),length(udelay));
        IL6tsevav = cell(length(uimax),length(udelay));
        Vtshedav = cell(length(uimax),length(udelay));
        diff_Vtshedav = cell(length(uimax),length(udelay));
        for ii = 1:length(udelay)
            for jj = 1:length(uimax)
                condti = (Delay(:) == udelay(ii)) & (Imax(:) == uimax(jj));
                bufV = log10(V_nomtav_merck(2,condti)) - log10(V_treatavt_merck(2,condti));
                diff_Vav{jj,ii} = bufV;

                bufIL6 = (IL6_nomtav_merck(2,condti)) - (IL6_treatavt_merck(2,condti));
                bufIL6 = bufIL6*mw.il6;
                bufIL6(bufIL6<=0.01) = 0.01; % lower bound of IL-6 differences
                diff_IL6av{jj,ii} = bufIL6;
                bufIL6 = (max(IL6_nomav_merck(:,condti)) - max(IL6_treatav_merck(:,condti)))*mw.il6;
                bufIL6(bufIL6<0.0) = NaN;
                bufIL6(bufIL6<=0.01) = 0.01; % lower bound of max(IL-6) differences
                diff_IL6maxav{jj,ii} = bufIL6;

                bufIL6 = max(IL6_treatav_merck(:,condti))*mw.il6;
                bufIL62 = (max(IL6_nomav_merck(:,condti)) - max(IL6_treatav_merck(:,condti)))*mw.il6;
                IL6maxav{jj,ii} = bufIL6;


                % NOTE - Time spent in servere condition not used in main results yet
                tIL6sevav = zeros(length(find(condti)),1);
                tvshedav = NaN(length(find(condti)),1);
                tvshednom = NaN(length(find(condti)),1);
                for kk = 1:length(find(condti))
                    condind = find(condti);
                    t1il6 = T_sampleav(find(IL6_treatav_merck(:,condind(kk))*mw.il6>70,1,'first'))/24;
                    
                    t1v =  tvmax_3s_merck(kk)+udelay(ii);           
                    t2v = T_sampleav(find(V_treatav_merck(:,condind(kk))<=1e2))/24;
                    t2v = t2v(find(t2v>t1v,1,'first'));
                    t2v = min(t2v,25);

                    t2vnom = T_sampleav(find(V_nom_merck(:,condind(kk))<=1e2))/24;
                    t2vnom = t2vnom(find(t2vnom>t1v,1,'first'));
                    t2vnom = min(t2vnom,25);
            
                    if ~isempty(t1il6)
                        t2il6 = T_sampleav(find(IL6_treatav_merck(:,condind(kk))*mw.il6>70,1,'last'))/24;
                        tIL6sevav(kk) = t2il6 - t1il6;
                    end

                    if ~isempty(t2v)
                        tvshedav(kk) = t2v - t1v;
                    end

                    if ~isempty(t2vnom)
                        tvshednom(kk) = t2vnom - t1v;
                    end
                end
                tIL6sevav(bufIL62<0.0) = NaN;
                tIL6sevav(tIL6sevav>30) = 30;
                IL6tsevav{jj,ii} = tIL6sevav;
                Vtshedav{jj,ii} = tvshedav;
                diff_Vtshedav{jj,ii} = (tvshednom - tvshedav)./tvshednom;
            end
        end


        % SECTION - Difference in viral load: average PBO - averge TRT
            % Mean with bootstrapped 99th percentile interval of mean
            pbo = V_nomtav_merck([4,6],:);
            pbo(pbo<=1) = 1; % in order to set log(V)<0 = 0
            trt = V_treatavt_merck([4,6],:);
            trt(trt<=1) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVav_bs = -1*ones(n_bs,size(pbo,1));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(samples_merck),length(samples_merck),true);
                pbo_bs = pbo(:,ind_bs);
                trt_bs = trt(:,ind_bs);
                diff_meanVav_bs(kk,:) = (mean(log10(pbo_bs),2) - mean(log10(trt_bs),2))';
            end

            diff_meanVav_merck.(mm{1}) = diff_meanVav_bs;
            mean(diff_meanVav_bs)

        % !SECTION - end Difference in viral load: average PBO - averge TRT

        % SECTION Calulate relative risk reduction (RRR) in severity end-point
            [maxil6,maxiil6] = max(squeeze(state_arraynom_merck(:,2,:))*mw.il6); % Maximum IL-6 levels in PBO group
            il6thresh = 40; % IL-6 threshold for events
            indsevere = find(maxil6>il6thresh);
            nsevere = length(find(maxil6>il6thresh));

            [maxil6_treat,maxiil6_treat] = max(squeeze(simav_merck.(mm{1})(:,2,:))*mw.il6);
            find(maxil6_treat>il6thresh)

            % Bootstrap IL-6 reduction after anti-viral 
            n_bs = 1000;
            rrr_av_bs_merck = -1*ones(1000,1);
            pbo_av_bs_merck = -1*ones(1000,1);
            trt_av_bs_merck = -1*ones(1000,1);
            for kk = 1:n_bs
                % NOTE - Bootsrapping to obtain C.I. around severity predictions.
                bsind = randsample(1:length(samples_merck),length(samples_merck),true);
                bsind2 = randsample(indsevere,nsevere,true); 
                event_pbo = find(maxil6(bsind)>=40);
                pbo_av_bs_merck(kk) = length(event_pbo)/length(samples_merck);
                trt_av_bs_merck(kk) = length(find(IL6maxav{ii,1}(bsind2)>il6thresh))/length(samples_merck); % Events in TRT arm
                rrr_av_bs_merck(kk) = 1 - length(find(IL6maxav{ii,1}(bsind2)>il6thresh))/nsevere; % RRR
            end
            pbo_av_trial_merck.(mm{1}) = pbo_av_bs_merck;
            trt_av_trial_merck.(mm{1}) = trt_av_bs_merck;
            rrr_av_trial_merck.(mm{1}) = rrr_av_bs_merck;
    end
    % !SECTION - end Calulate relative risk reduction (RRR) in severity end-point

% !SECTION - end Analyze solutions/Compare PBO and TRT


%% SECTION - Generate Plots for anti-viral molnupiravir simulations

    % SECTION - Mean change from baseline viral load & mean viral load
        mean_diff_nombase = mean(log10(Vtime.x800mg.V_nomtav_merck) - log10(Vtime.x800mg.V_nomtav_merck(1,:)),2);
        mean_diff_treatbase = mean(log10(Vtime.x800mg.V_treatavt_merck) - log10(Vtime.x800mg.V_treatavt_merck(1,:)),2);
        std_diff_nombase = std(log10(Vtime.x800mg.V_nomtav_merck) - log10(Vtime.x800mg.V_nomtav_merck(1,:)),[],2);
        std_diff_treatbase = std(log10(Vtime.x800mg.V_treatavt_merck) - log10(Vtime.x800mg.V_treatavt_merck(1,:)),[],2);
        figure,
        tiledlayout(2,2),
        nexttile,
        pl1 = errorbar([0,2,3,4,5,6,10]+1,mean_diff_nombase,std_diff_nombase/sqrt(62),'LineWidth',2,'Color','k');
        hold on, 
        pl2 = errorbar([0,2,3,4,5,6,10]+1,mean_diff_treatbase,std_diff_treatbase/sqrt(62),'LineWidth',2);
        xlim([1 8])
        ylim([-5 0])
        ylabel('Mean Change from Baseline Viral Load log_{10} copies/mL')
        xlabel('Days from start of treatment')
        hold on, 
        plot(molnupiravir_dbsl.Time(1:end-1), molnupiravir_dbsl.Placebo(1:end-1),'Color',[0,0,0],'LineStyle','--','LineWidth',2,'Marker','+','MarkerSize',20)
        plot(molnupiravir_dbsl.Time(1:end-1), molnupiravir_dbsl.x800mg(1:end-1),'Color',[0,0.4470,0.7410],'LineStyle','--','LineWidth',2,'Marker','+','MarkerSize',20)
        ax = gca;
        ax.FontSize = 14;
        grid on
        lg = legend('Virtual population PBO (mean \pm S.E.)','Virtual population Treated (mean \pm S.E.)','Merck Ph II Trial PBO','Merck Ph II Trial Treated','Location','SouthWest');
        lg.Title.String = '4d post symptom onset';

        nexttile,
        errorbar([0,2,3,4,5,6,10]+1,mean(log10(Vtime.x800mg.V_nomtav_merck),2),std(log10(Vtime.x800mg.V_nomtav_merck),[],2)/sqrt(62),'LineWidth',2)
        hold on, 
        errorbar([0,2,3,4,5,6,10]+1,mean(log10(Vtime.x800mg.V_treatavt_merck),2),std(log10(Vtime.x800mg.V_treatavt_merck),[],2)/sqrt(62),'LineWidth',2)
        xlim([1 8])
        ylim([0 7])
        ylabel('Mean Viral Load log_{10} copies/mL')
        xlabel('Day of trial')
        hold on, 
        scatter(molnupiravir_vl.Time, molnupiravir_vl.Placebo,50,'filled')
        scatter(molnupiravir_vl.Time, molnupiravir_vl.x800mg,50,'filled')
        ax = gca;
        ax.FontSize = 14;
        grid on
        lg = legend('Simulated PBO','Simulated Treated','Trial PBO','Trial Treated','Location','SouthWest');
        lg.Title.String = '4d post symptom onset';


        nexttile,
        barp = bar([[0.547;0.534],[mean(diff_meanVav_merck.x800mg(:,1));mean(diff_meanVav_merck.x800mg(:,2))]]);
        hold on,
        e2 = errorbar(barp(2).XEndPoints(2),mean(diff_meanVav_merck.x800mg(:,2)), (mean(diff_meanVav_merck.x800mg(:,2)) - prctile(diff_meanVav_merck.x800mg(:,2),0.0)),(prctile(diff_meanVav_merck.x800mg(:,2),100) - mean(diff_meanVav_merck.x800mg(:,2))),'LineWidth',1.5,'LineStyle','none','Color','k');   
        e2 = errorbar(barp(2).XEndPoints(1),mean(diff_meanVav_merck.x800mg(:,1)), (mean(diff_meanVav_merck.x800mg(:,1)) - prctile(diff_meanVav_merck.x800mg(:,1),0.0)),(prctile(diff_meanVav_merck.x800mg(:,1),100) - mean(diff_meanVav_merck.x800mg(:,1))),'LineWidth',1.5,'LineStyle','none','Color','k');   
        e1 = errorbar(barp(1).XEndPoints(2),0.534, (0.534 - 0.157), (0.91 - 0.534),'LineWidth',1.5,'LineStyle','none','Color','k');   
        e1 = errorbar(barp(1).XEndPoints(1),0.547, (0.547 - 0.159), (0.935 - 0.547),'LineWidth',1.5,'LineStyle','none','Color','k');   
        grid on
        xlabel('[Treatment starts at Day 1]')
        ylabel('log_{10}[reduction in viral load]')
        xticklabels({'Day 5','Day 7'})
        ylim([0 2.0])
        set(gca,'FontSize',16)
        grid on
        leg = legend([barp,e1],'Molnupiravir 800mg Ph II','Model simulation');
        set(gcf, 'Position',  [300   527   1000   800])

    %!SECTION - end - Mean change from baseline viral load

% !SECTION - end Generate Plots for anti-viral molnupiravir simulations