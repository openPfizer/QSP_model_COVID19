% script to generate all the figures from manuscript

addpath(genpath('.'))
run plausible_figure.m
run BE_Blaze1Ph3.m
run regen_cov.m
run severity_plot.m
run dVL_RRR_timing.m
run supplementary_figures.m

h =  findobj(allchild(0),'flat','Type','figure');
for ii = 1:length(h)
    figno = h(ii);
    figname2 = get(figno,'Name');
    print(figno,figname2,'-dpng');
    
end
close all

exit
