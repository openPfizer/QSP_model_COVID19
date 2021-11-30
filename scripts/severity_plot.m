%% NOTE: Script to plot severity improvements for both Blaze-1 Ph3 trial and REGEN-COV Ph3 trial

%% SECTION: Plot event rates
    figure,
        scatter([1,3,8,10],[7,2.1,4.6,1.3],512,repmat([0,0,0; 0,0.4470,0.7410],2,1),'filled','Marker','d')
        hold on,
        scat1 = scatter(0,-5,512,'filled','Marker','d');
        scat1.CData = [0,0,0];
        scat2 = scatter(0,-5,512,'filled','Marker','d');
        scat2.CData = [0,0.4470,0.7410];
        vsplot = violinplot([-1;100*(pbo_ly);-1;100*(trt_ly);-1;-1;-1;-1;100*(pbo_regn);-1;100*(trt_regn)], [0;1*ones(size(pbo_ly));2;3*ones(size(trt_ly));4;5;6;7;8*ones(size(pbo_regn));9;10*ones(size(pbo_ly))],'ViolinAlpha',1,'ViolinColor',[0.8500,0.3250,0.0980],'BoxColor',[0,0,0],'ShowData',false);
        vsplot(2).ViolinColor = [0 0 0];
        vsplot(2).ViolinAlpha = 0.4;
        vsplot(4).ViolinPlot.FaceColor = [0.8500,0.3250,0.0980];
        vsplot(4).ViolinColor = [0.8500,0.3250,0.0980];
        vsplot(9).ViolinPlot.FaceColor = [0 0 0];
        vsplot(9).ViolinColor = vsplot(2).ViolinColor;
        vsplot(9).ViolinAlpha = vsplot(2).ViolinAlpha;
        vsplot(11).ViolinPlot.FaceColor = [0.8500,0.3250,0.0980];
        vsplot(11).ViolinColor = [0.8500,0.3250,0.0980];
        ylim([0 11])
        xlim([-1 12])
        xticks([2.5, 9.5])
        xticklabels({'Blaze-1 Ph3','REGEN-COV Ph3'})
        row1 = {'Blaze-1 Ph3 [2800mg bamlanivimab','REGEN-COV Ph3'};
        row2 = {' + 2800mg etesevimab]', '[2.4g REGEN-COV]'};
        labelArray = [row1; row2]; 
        labelArray = strjust(pad(labelArray),'center');
        tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
        ax = gca(); 
        ax.XTickLabel = tickLabels;
        ax.FontSize = 14;
        ylabel('Event Rate (%)')
        legend([scat1,scat2,vsplot(2).ViolinPlot,vsplot(4).ViolinPlot],'Observed placebo','Observed treated group','Simulated placebo','Simulated treatment','Location','NorthEast')
        grid on
% !SECTION - end Plot event rates

%% SECTION: Plot RRR
    figure,
        plt = scatter([1,7],[70,69],512,[0,0.4470,0.7410],'filled','Marker','d');
        hold on
        ylabel('Relative risk reduction (%)')
        hold on,
        vsplot = violinplot([-1;100*rrr_ly;-1;-1;-1;-1;-1;100*rrr_regn],[0;1*ones(size(rrr_ly));2;3;4;5;6;7*ones(size(rrr_regn))],'ShowData',false,'ViolinAlpha',1);
        vsplot(2).ViolinColor = [0.8500,0.3250,0.0980];
        vsplot(2).ViolinPlot.FaceColor = [0.8500,0.3250,0.0980];
        vsplot(8).ViolinColor = [0.8500,0.3250,0.0980];
        vsplot(8).ViolinPlot.FaceColor = [0.8500,0.3250,0.0980];
        e2 = errorbar(7,69,69-46.2,82.1-69,'LineWidth',1.5,'Color','k');
        grid on
        xticks([1.5 7.5])
        row1 = {'Blaze-1 Ph3 [2800mg bamlanivimab','REGEN-COV Ph3'};
        row2 = {' + 2800mg etesevimab]', '[2.4g REGEN-COV]'};
        labelArray = [row1; row2]; 
        labelArray = strjust(pad(labelArray),'center');
        tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
        ax = gca(); 
        ax.XTickLabel = tickLabels;
        ax.FontSize = 14;
        ax.FontSize = 14;
        ylim([0 100])
        xlim([-1 9])
        legend([plt,vsplot(2).ViolinPlot],'Observed','Simulated','Location','Best')
% !SECTION - end Plot RRR