%% Merge data for fitting
function [merged_data, statetbl, varargout] = merge_optimdata(data,data_dictionary,varargin)

    task = varargin{1};
    switch task
        case "merge"

            values = [data_dictionary.species_names(1:end-4,2);"CD4";'Monocytes'];
            keys = ["V"
                    "AT1"
                    "AT2"
                    "I"
                    "dAT1"
                    "dAT2"
                    "l_pDC"
                    "l_M1"
                    "l_N"
                    "l_Th1"
                    "l_Th17"
                    "l_CTL"
                    "l_Treg"
                    "l_TNFa"
                    "l_IL6"
                    "l_IL1b"
                    "l_IFNb"
                    "l_IFNg"
                    "l_IL2"
                    "l_IL12"
                    "l_IL17"
                    "l_IL10"
                    "l_TGFb"
                    "l_GMCSF"
                    "l_SPD"
                    "l_FER"
                    "TNFa"
                    "IL6"
                    "IL1b"
                    "IFNa2"
                    "IFNg"
                    "IL2"
                    "IL12"
                    "IL17"
                    "IL10"
                    "TGFb"
                    "GMCSF"
                    "SPD"
                    "FER"
                    "Ab"
                    "pDC"
                    "M1"
                    "N"
                    "Th1"
                    "Th17"
                    "Activated_CD8"
                    "Treg"
                    "CRPExtracellular"
                    "CRP"
                    "pDC1"
                    "Activated_CD4"
                    "Monocytes"];

            match_st = containers.Map(keys, values);
            st_match = containers.Map(values, keys);

            if varargin{2} == "states"
                optim_states = varargin{3};
            else
                optim_states = ["V"
                            "TNFa_c"
                            "IL6_c"
                            "IL1b_c"
                            "IFNb_c"
                            "IFNg_c"
                            "IL2_c"
                            "IL12_c"
                            "IL10_c"
                            "GMCSF_c"
                            "N_c"
                            "pDC_c"
                            "CTL_c"
                            "Blood_CRP"
                            "CD4"];
            end

                maxbound.V = Inf;
                maxbound.TNFa_c = 1000;
                maxbound.IL6_c = 2000;
                maxbound.IFNg_c = 500;
                maxbound.IFNb_c = 5000;
                maxbound.IL1b_c = 1000;
                maxbound.IL2_c = 500;
                maxbound.IL12_c = 500;
                maxbound.IL10_c = 800;
                maxbound.GMCSF_c = 1600;
                maxbound.pDC_c = 30;
                maxbound.N_c = 5000;
                maxbound.CD4 = 50;
                maxbound.CTL_c = 50;
                maxbound.Blood_CRP = 300;




            optim_dataset = {["gastine"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["lucas"]
                            ["mudd"]
                            ["lucas"]
                            ["lucas"]
                            ["mudd"]
                            ["manson"]
                            ["mudd"]
                            };

            func = @(x) [median(x, 'omitnan'), mean(x, 'omitnan'), max(x, [],'omitnan'), min(x, [],'omitnan')];

            kk = 1;
            for ii = optim_states'

                bufstate = [];
                buftime = [];
                bufcond = [];
                for jj = optim_dataset{kk}
                    state = st_match(ii);
                        if strcmp(state,"V")
                            data.(jj).Condition(data.(jj).Condition == "Moderate" | data.(jj).Condition == "Severe") = "Select";
                            cond1v = data.(jj).vl_quality == 1;
                            cond2v = data.(jj).sample_site == "nasopharyngeal";
                            cond3v = data.(jj).drug_quality == 1;
                            % data.(jj) = table(data.Time(cond1v & cond2v & cond3v & cond4v), data.(state), data.Condition,'VariableNames',{'Time','Condition',match_st(state)});
                            cond_opt = (cond1v & cond2v & cond3v);
                        else
                            data.(jj).Condition(data.(jj).Condition == "Severe"  |  data.(jj).Condition == "Moderate"  |data.(jj).Condition == "COVID-19") = "Select";
                            cond_opt = ones(size(data.(jj).Condition));
                        end
                        cond1 = data.(jj).(state)<=maxbound.(ii);
                        cond2 = data.(jj).Time<20;
                        cond3 = data.(jj).Condition == "Select" | data.(jj).Condition == "Healthy";
                        bufstate = [bufstate;data.(jj).(state)(cond1 & cond2 & cond3 & cond_opt)];
                        buftime = [buftime; data.(jj).Time(cond1 & cond2 & cond3 & cond_opt)];
                        bufcond = [bufcond; data.(jj).Condition(cond1 & cond2 & cond3 & cond_opt)];
                end
                kk = kk+1;
                tbl = table(buftime,bufcond,bufstate,'VariableNames',{'Time','Condition',match_st(state)});
                statetbl.(ii) = tbl;
                dat.(ii) = varfun(func,tbl,'GroupingVariables',{'Time','Condition'},'InputVariables',(ii));
                dat.(ii).Properties.VariableNames(4) = "stat";
                dat.(ii) = splitvars(dat.(ii),'stat','NewVariableNames',{'median','mean','max','min'});
                    
            end

            merged_data = dat;
        
        case "plot"
            if strcmp(varargin{1},"plot")

                datat0 = varargin{2};
                datap0 = varargin{3};
                datap.mean = -1*ones(length(unique(datat0)),1);
                datap.min = -1*ones(length(unique(datat0)),1);
                datap.max = -1*ones(length(unique(datat0)),1);
                datat = unique(datat0);
                for ii = 1:length(unique(datat0))
                    datap.mean(ii) = mean(datap0(datat0 == datat(ii)),'omitnan');
                    datap.min(ii) = min(datap0(datat0 == datat(ii)));
                    datap.max(ii) = max(datap0(datat0 == datat(ii)));
                end

                varargout{1} = datat;
                varargout{2} = datap;
                merged_data = []; 
                statetbl = [];
            
            end

    end


end

% cond1 = data.lucas.Time == 23;
% cond2 = data.lucas.Condition == "Moderate";

% data.lucas.V(cond1)
% fnames = fieldnaems(mrgd)



%     for ii = fieldnames(mrgd)'
%         kk = [ismember(ii,"V")
%               ismember(ii,cen_cytokines)
%               ismember(ii,cen_cells)];

%         switch jj(kk)
            
%             writedata = mrgd.(ii{:});
%             case "V"
%                 writedata.dataunits = strings(length(writedata.Time),1)+,"cp/mL";
%                 writedata.modelunits = strings(length(writedata.Time),1)+"cp/mL";
%             case "cen_cytokines"
%                 writedata.dataunits = strings(length(writedata.Time),1)+"pg/mL";
%                 writedata.modelunits = strings(length(writedata.Time),1)+"pmol/mL";
%             case "cen_cells"
%                 writedata.dataunits = strings(length(writedata.Time),1)+"cells/uL";
%                 writedata.modelunits = strings(length(writedata.Time),1)+"cells/mL";
%         end

%         
%     end

% state = "IL6_c";

% % cond1 = mrgd.(state).Condition == "COVID-19";
% % cond2 = mrgd.(state).Condition == "Moderate"
% figure, 
% plot(mrgd.(state).Time, mrgd.(state).stat(:,1))
% hold on, 
% plot(mrgd.(state).Time, mrgd.(state).stat(:,3))
% plot(mrgd.(state).Time, mrgd.(state).stat(:,4))
% set(gca, 'YScale', 'log')

% set(gca, 'YScale', 'log')
% ylabel(state)
% xlabel('Days')
% legend("Severe","Moderate")

% cond1 = mrgd.(state).Condition == "Severe";
% cond2 = mrgd.(state).Condition == "Moderate";
% figure, 
% scatter(mrgd.(state).Time(cond1), mrgd.(state).stat(cond2,2))
% hold on, 
% scatter(mrgd.(state).Time(cond2), mrgd.(state).stat(cond1,2))
% set(gca, 'YScale', 'log')
% ylabel(state)
% xlabel('Days')
% legend("Severe","Moderate")
% prctile(data.lucas.IFNg, 95)
% length(find(data.lucas.IFNg>500))

% cond1 = data.lucas.Condition == "Severe";
% plot_data('data',data,'cytokine',' ','cell',"pDC",'studies',"lucas")
% set(gca, 'YScale','log')
% plot_data('data',data,'cytokine','IFNg','studies',"mann")
% set(gca, 'YScale','log')
% prctile(data.lucas.TNFa(cond1), 90)
