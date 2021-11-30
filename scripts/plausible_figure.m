% NOTE: Script to plot plausible population against cytokine and viral load data

load("plausible_fig.mat")

%% SECTION - Plot plausible population wrt inferred time of infection (assumed incubation period in data = 4.5d)

        data_dictionary_orig = get_data_dictionary(); % load model dictionary
        mw = data_dictionary_orig.mw;
        
        figure,
        subplot(5,2,3)
        ax = gca;
        state = "IFNb_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.ifnb*16/20;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')
        
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'Type I IFN (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on,

        subplot(5,2,2)
        ax = gca;
        state = "IL6_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.il6;
        statehigh = prctile(statep,100,2);
        statelow = prctile(statep,0.0,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)
    
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-6 (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,4)
        ax = gca;
        state = "IFNg_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.ifng;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)   

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IFNg (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,5)
        ax = gca;
        state = "TNFa_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.tnfa;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)   
       
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'TNFa (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,6)
        ax = gca;
        state = "IL10_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.il10;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)  

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-10 (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,7)
        ax = gca;
        state = "IL1b_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.il1b;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)   

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-1b (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,8)
        ax = gca;
        state = "IL12_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.il12;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2)   
    
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-12 (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,9)
        ax = gca;
        state = "IL2_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:))*mw.il2;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)% set(gca,'YScale','log')
        hold on
        plot(T/24,mean(statep,2),'LineWidth',2) 

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).healthy.infectionTh,pdata.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata.(state).covid.infectionTd,pdata.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-2 (pg/mL)'})
        xlabel('Days from infection')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on
        lgd1 = legend('99^{th} percentile interval (Simulated)','Mean (Simulated)','Mean \pm range clinical data (Healthy)','Mean \pm range clinical data (COVID-19)','Position',[0.1965+0.4 0.1495 0.2600 0.0784]);

        subplot(5,2,1)
        state = "V";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = squeeze(state_arrayplaus1(:,iplot,:));
        statep(statep<1) = 1;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        patch([T/24;flipud(T/24)],log10([statehigh;flipud(statelow)])', 'k','FaceAlpha',0.2)% set(gca,'YScale','log')
        hold on
        plot(T/24,mean(log10(statep),2),'LineWidth',2) 

        [mrgd1,~] = merge_optimdata(pdata,data_dictionary_orig,'merge','states',"V");
        errorbar(mrgd1.(state).Time, (mrgd1.(state).mean),(mrgd1.(state).mean-mrgd1.(state).min),(-mrgd1.(state).mean+mrgd1.(state).max),'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        xlabel('Days from infection')
        ylabel ({'log_{10} Viral', 'RNA copies/mL'})
        set(gca,'FontSize',14)
        set(gcf,'Position',  [150   527   1000   1000])
        grid on

% !SECTION - end Plot plausible population wrt inferred time of infection (assumed incubation period in data = 4.5d)

%% SECTION - Plot plausible population wrt assumed incubation period for each vsub (time of simulated symptom onset = time of simulated viral load peak)

        T = T_sample;
        Tind = find(T<20*24,1,'last');
        T = T(freq);
        figure,
        subplot(5,2,3)
        ax = gca;
        state = "IFNb_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.ifnb*16/20, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)
        

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')
        
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'Type I IFN (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on,

        subplot(5,2,2)
        ax = gca;
        state = "IL6_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.il6, tvmax_3s_plaus, T);
        statehigh = prctile(statep,100,2);
        statelow = prctile(statep,0.0,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)% set(gca,'YScale','log')
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)
    
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-6 (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,4)
        ax = gca;
        state = "IFNg_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.ifng, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)% set(gca,'YScale','log')
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IFNg (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,5)
        ax = gca;
        state = "TNFa_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.tnfa, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)% set(gca,'YScale','log')
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)
       
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'TNFa (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,6)
        ax = gca;
        state = "IL10_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.il10, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2) 

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-10 (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,7)
        ax = gca;
        state = "IL1b_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.il1b, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-1b (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,8)
        ax = gca;
        state = "IL12_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.il12, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2)  
    
        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-12 (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on

        subplot(5,2,9)
        ax = gca;
        state = "IL2_c";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:))*mw.il2, tvmax_3s_plaus, T);
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],[statehigh;flipud(statelow)]', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(statep,2),'LineWidth',2) 

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).healthy.infectionTh,pdata2.(state).healthy.healthy);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','b','LineStyle','none','Marker','o')

        [~,~,datat, datap] = merge_optimdata([], data_dictionary_orig, 'plot',pdata2.(state).covid.infectionTd,pdata2.(state).covid.covid);
        errorbar(datat, datap.mean,datap.mean-datap.min,-datap.mean+datap.max,'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        ylabel({'IL-2 (pg/mL)'})
        xlabel('Days from symptom onset')
        hold on
        xlim([0 20])
        ticks = 10.^(floor(log10(min(statelow))):ceil(log10(max(statehigh))));
        yticks(ticks);
        set(ax, 'YScale', 'log')
        set(gca,'FontSize',14)
        grid on
        lgd2 = legend('99^{th} percentile interval (Simulated)','Mean (Simulated)','Mean \pm range clinical data (Healthy)','Mean \pm range clinical data (COVID-19)','Position',[0.1965+0.4 0.1495 0.2600 0.0784]);

        subplot(5,2,1)
        ax = gca;
        state = "V";
        [~,iplot] = intersect(statenames_to_store,state);
        statep = adjust_tfso(squeeze(state_arrayplaus1(:,iplot,:)), tvmax_3s_plaus, T);
        statep(statep<1) = 1;
        statehigh = prctile(statep,99.5,2);
        statelow = prctile(statep,0.5,2);
        Tplot = [0,2,4,6,7,9,10,12,15,20]';
        patch([Tplot;flipud(Tplot)],log10([statehigh;flipud(statelow)])', 'k','FaceAlpha',0.2)
        hold on
        plot(Tplot,mean(log10(statep),2),'LineWidth',2) 

        [mrgd1,~] = merge_optimdata(pdata2,data_dictionary_orig,'merge','states',"V");
        errorbar(mrgd1.(state).Time, (mrgd1.(state).mean),(mrgd1.(state).mean-mrgd1.(state).min),(-mrgd1.(state).mean+mrgd1.(state).max),'LineWidth',2,'Color','k','LineStyle','none','Marker','o')
        xlabel('Days from symptom onset')
        ylabel ({'log_{10} Viral', 'RNA copies/mL'})
        set(gca,'FontSize',14)
        set(gcf,'Position',  [150   527   1000   1000])
        grid on

% !SECTION - Plot plausible population wrt assumed incubation period for each vsub (time of simulated symptom onset = time of simulated viral load peak)