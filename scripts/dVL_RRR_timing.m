%% NOTE: Plot sensitivity of reduction in viral and RRR to time of intervention

%% SECTION - Run time of intervention simulations
    udelay = [1, 2, 4, 6, 8, 10, 12]'; % administration time with respect to exposure
    filter_3s = samples_indmin;
    delay = repmat(repelem(udelay,1,size(filter_3s,1)),1,1); % array of delays for each sample
    Delay = delay';
    n_vpav = numel(Delay);
    n_freqav = 12;
    n_time_pointsav = length(0:0.1:(40*24)/n_freqav); % tstart : tstep: tstop divided by sampling frequency
    state_arraybnab_time = zeros(n_time_pointsav,2,numel(Delay));
    tmp_parameters = data_dictionary_orig.parameters;
    err_vectorbnab = -1*ones(numel(Delay),1);
    parfor ii = 1:n_vpav
        sample_indexn = ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1));
        % updating model dictionary with specific parameter sample
        [data_dictionary, virus_innoculation] = update_parameters_ext(data_dictionary_orig,tmp_parameters,sample_indexn,perturbation_name_vector,samples_indmin); % Update data_dictionary for each sample
        data_dictionary.parameters.tau_Ab = 24*(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii));
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
        state_arraybnab_time(:,:, ii) = X(1:n_freqav:end,[1,28]);


        if size_sample < n_time_pointsav
            err_vectorbnab = [err_vectorbnab;ii];
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

% !SECTION - end Run time of intervention simulations

%% SECTION - Analyze solutions/Compare PBO and TRT for nAbs 
        mw = data_dictionary_orig.mw; % molecular weights for species
        % state_arraynomav = state_arraynom;
        V_nomab = squeeze(state_arraynomav(:,1,:)); % get viral load
        V_nomab = repmat(V_nomab,1,numel(Delay)/size(filter_3s,1)); % adjust to same dimensions as nAb solutions
        V_treatab = squeeze(state_arraybnab_time(:,1,:)); % treated viral load
        V_nomtab = -1*ones(8,numel(Delay)); % initialize array to store solutions at trial time points
        V_treatabt = -1*ones(8,numel(Delay));

        % Initialize IL-6 arrays
        IL6_nomab = squeeze(state_arraynomav(:,2,:));
        IL6_nomab = repmat(IL6_nomab,1,numel(Delay)/size(filter_3s,1));
        IL6_treatab = squeeze(state_arraybnab_time(:,2,:));
        IL6_nomtav = -1*ones(8,numel(Delay));
        IL6_treatabt = -1*ones(8,numel(Delay));
        IL6_nomabmax = -1*ones(1,numel(Delay));
        % Get viral load and IL-6 at trial time points
        for ii = 1:numel(Delay)
            % Get time vector indices for trial time points
            ind = [find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii)),1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+2,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+4,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+5,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+6,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+7,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+10,1,'last');
            find(T_sampleav/24<=(tvmax_3s(ii-size(filter_3s,1)*floor((ii-1)/size(filter_3s,1)))+Delay(ii))+15,1,'last')];
            V_nomb = V_nomab(ind,ii);
            V_nomb(V_nomb<=1e0) = 1; % in order to set log(V)<0 = 0
            V_nomtab(:,ii) = V_nomb;

            V_treatabb = V_treatab(ind,ii);
            V_treatabb(V_treatabb<=1e0) = 1; % in order to set log(V)<0 = 0
            V_treatabt(:,ii) = V_treatabb;

            IL6_nomb = IL6_nomab(ind,ii);
            IL6_nomtav(:,ii) = IL6_nomb;

            IL6_treatabb = IL6_treatab(ind,ii);
            IL6_treatabt(:,ii) = IL6_treatabb;
            
        end

        % Calculate differences in viral load and IL-6 for PBO_i vs TRT_i (paired)
        udelay = unique(delay);
        diff_Vab = cell(1,length(udelay));
        diff_IL6ab = cell(1,length(udelay));
        diff_IL6maxab = cell(1,length(udelay));
        IL6maxab = cell(1,length(udelay));
        IL6tsevab = cell(1,length(udelay));
        Vtshedab = cell(1,length(udelay));
        diff_Vtshedab = cell(1,length(udelay));
        for ii = 1:length(udelay)
                condti = (Delay(:) == udelay(ii));
                bufV = log10(V_nomtab(2,condti)) - log10(V_treatabt(2,condti));
                diff_Vab{1,ii} = bufV;

                bufIL6 = (IL6_nomtav(2,condti)) - (IL6_treatabt(2,condti));
                bufIL6 = bufIL6*mw.il6;
                bufIL6(bufIL6<=0.01) = 0.01; % lower bound of IL-6 differences
                diff_IL6ab{1,ii} = bufIL6;
                bufIL6 = (max(IL6_nomab(:,condti)) - max(IL6_treatab(:,condti)))*mw.il6;
                bufIL6(bufIL6<0.0) = NaN;
                bufIL6(bufIL6<=0.01) = 0.01; % lower bound of max(IL-6) differences
                diff_IL6maxab{1,ii} = bufIL6;

                bufIL6 = max(IL6_treatab(:,condti))*mw.il6;
                bufIL62 = (max(IL6_nomab(:,condti)) - max(IL6_treatab(:,condti)))*mw.il6;
                IL6maxab{1,ii} = bufIL6;

                % Calculate time spent in servere condition (IL-6>70 pg/mL); max time allowed in severe condition before "mortality" = 25 days. min(t2v,25).
                % NOTE - Time spent in servere condition not used in main results yet
                tIL6sevav = zeros(length(find(condti)),1);
                tvshedav = NaN(length(find(condti)),1);
                tvshednom = NaN(length(find(condti)),1);
                for kk = 1:length(find(condti))
                    condind = find(condti);
                    t1il6 = T_sampleav(find(IL6_treatab(:,condind(kk))*mw.il6>70,1,'first'))/24;
                    
                    t1v =  tvmax_3s(kk)+udelay(ii);           
                    t2v = T_sampleav(find(V_treatab(:,condind(kk))<=1e2))/24;
                    t2v = t2v(find(t2v>t1v,1,'first'));
                    t2v = min(t2v,25);

                    t2vnom = T_sampleav(find(V_nomab(:,condind(kk))<=1e2))/24;
                    t2vnom = t2vnom(find(t2vnom>t1v,1,'first'));
                    t2vnom = min(t2vnom,25);
            
                    if ~isempty(t1il6)
                        t2il6 = T_sampleav(find(IL6_treatab(:,condind(kk))*mw.il6>70,1,'last'))/24;
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
                IL6tsevab{1,ii} = tIL6sevav;
                Vtshedab{1,ii} = tvshedav;
                diff_Vtshedab{1,ii} = (tvshednom - tvshedav)./tvshednom;
        end


        % SECTION - Difference in viral load: average PBO - averge TRT
            % Mean with bootstrapped 99th percentile interval of mean
            pbo = reshape(V_nomtab(5,:),length(samples_indmin),length(udelay));
            pbo(pbo<=1e0) = 1; % in order to set log(V)<0 = 0
            trt = reshape(V_treatabt(5,:),length(samples_indmin),length(udelay));
            trt(trt<=1e0) = 1; % in order to set log(V)<0 = 0
            n_bs = 1000;
            diff_meanVab_bs = -1*ones(n_bs,length(udelay));
            for kk = 1:n_bs
                ind_bs = randsample(1:length(samples_indmin),length(samples_indmin),true);
                pbo_bs = pbo(ind_bs,:);
                trt_bs = trt(ind_bs,:);
                diff_meanVab_bs(kk,:) = (mean(log10(pbo_bs),1) - mean(log10(trt_bs),1))';
            end
            
        % !SECTION - end Difference in viral load: average PBO - averge TRT


        % SECTION Calulate relative risk reduction (RRR) in severity end-point
            [maxil6,maxiil6] = max(squeeze(state_arraynomav(:,2,:))*mw.il6); % Maximum IL-6 levels in PBO group
            il6thresh = 40; % IL-6 threshold for events
            indsevere = find(maxil6>il6thresh);
            nsevere = length(find(maxil6>il6thresh));

            % Bootstrap IL-6 reduction after anti-viral
            n_bs = 1000;
            rrr_nab = -1*ones(1000,length(udelay));
            pbo_nab = -1*ones(1000,1);
            trt_nab = -1*ones(1000,length(udelay));
            for kk = 1:n_bs
                % NOTE Bootsrapping to obtain C.I. around severity predictions
                bsind = randsample(1:length(samples_indmin),length(samples_indmin),true);
                bsind2 = randsample(indsevere,nsevere,true); 
                event_pbo = find(maxil6(bsind)>=40);
                pbo_nab(kk) = length(event_pbo)/length(samples_indmin);
                trt_nab(kk,:) = [length(find(IL6maxab{1,1}(bsind2)>il6thresh));length(find(IL6maxab{1,2}(bsind2)>il6thresh));length(find(IL6maxab{1,3}(bsind2)>il6thresh));length(find(IL6maxab{1,4}(bsind2)>il6thresh));length(find(IL6maxab{1,5}(bsind2)>il6thresh));length(find(IL6maxab{1,6}(bsind2)>il6thresh));length(find(IL6maxab{1,7}(bsind2)>il6thresh))]/length(samples_indmin); % Events in TRT arm
                rrr_nab(kk,:) = 1 - [length(find(IL6maxab{1,1}(bsind2)>il6thresh));length(find(IL6maxab{1,2}(bsind2)>il6thresh));length(find(IL6maxab{1,3}(bsind2)>il6thresh));length(find(IL6maxab{1,4}(bsind2)>il6thresh));length(find(IL6maxab{1,5}(bsind2)>il6thresh));length(find(IL6maxab{1,6}(bsind2)>il6thresh));length(find(IL6maxab{1,7}(bsind2)>il6thresh))]/nsevere; % RRR
            end

            mean(rrr_nab)
            max(rrr_nab)
            min(rrr_nab)

            mean(diff_meanVab_bs)

        % !SECTION - end Calulate relative risk reduction (RRR) in severity end-point

% !SECTION - end Analyze solutions/Compare PBO and TRT for nAbs

%% SECTION - Generate figures for time of intervention of nAb 
figure, 
    subplot(1,2,1)
    errorbar(udelay, mean(diff_meanVab_bs), mean(diff_meanVab_bs) - prctile(diff_meanVab_bs,0.5),-mean(diff_meanVab_bs) + prctile(diff_meanVab_bs,99.5), 'LineWidth',2,'Color','k')
    ylabel('log_{10}[reduction in viral load]')
    xlabel('Intervention relative to symptom onset (days)')
    xlim([1 13])
    ylim([0.1 1.6])
    grid on
    set(gca,'FontSize',14)

    subplot(1,2,2)
    errorbar(udelay, mean(rrr_nab)*100, (mean(rrr_nab) - prctile(rrr_nab,0))*100,(-mean(rrr_nab) + prctile(rrr_nab,100))*100, 'LineWidth',2,'Color','k')
    ylabel({'Relative risk reduction'; 'in simulated events (%)'})
    xlabel('Intervention relative to time from symptom onset (days)')
    ylim([0 100])
    xlim([1 13])
    grid on
    set(gca,'FontSize',14)

    set(gcf, 'Position',  [150   527   1000   400])

% !SECTION - end Generate figures for time of intervention of nAb   