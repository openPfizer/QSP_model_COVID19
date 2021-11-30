%% script to generate all the figures from manuscript
addpath(genpath('.'))
run potency_eua_v2.m

h =  findobj(allchild(0),'flat','Type','figure');
for ii = 1:length(h)
    figno = h(ii);
    figname2 = get(figno,'Name');
    figname2 = strcat("/covid_sequel/docs/",figname2)
    print(figno,figname2,'-dpng');
    
end
close all

exit