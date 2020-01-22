function [] = stability_plot(C,Time,S,N,VI)
% plot Markov stability curve
figure_handle = figure;
set(0,'CurrentFigure',figure_handle);
t = length(Time);
subplot(4,1,1), ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
xlabel('Markov time');
set(ax(1),'YScale','log');
set(ax(2),'YScale','log');
set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(ax(2),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(get(ax(1),'Ylabel'),'String','Number of communities');
set(get(ax(2),'Ylabel'),'String','Stability');
set(ax(1),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [1 10^ceil(log10(max(N)))], 'XScale','log','XMinorGrid','on');
set(ax(2),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [10^floor(log10(min(S(N>1)))), 1], 'XScale','log');
ylabel('Number of communities');

subplot(4,1,2), semilogx(Time(1:t),VI(1:t),'Color',[0 0.7 0],'LineWidth',1.5);
set(gca, 'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], ...
    'YMinorGrid','on','XMinorGrid','on','YColor', [0 0.7 0]);
if max(VI)>0
    set(gca,'YLim', [0 max(VI)*1.1]);
end
xlabel('Markov time');
ylabel('Variation of information','Color',[0 0.7 0]);


[~,vi_mat] = varinfo(C',1);
%stability_plot(Time,t,S,N,VI,ComputeVI,figure_handle);
set(0,'CurrentFigure',figure_handle);
%subplot(4,1,[3,4])
% plot the vi(t,t')
ax1 = axes();
n = length(Time);
imagesc(ax1,[0,n-1],[0,n-1],vi_mat); colormap pink
axis(ax1,'off');
c = colorbar;
c.Location = 'southoutside';
c.Position = [0.66,0.066,0.2,0.02];
c.Label.String = 'VI(t,t'')';
% plot the Vi and N
ax2 = axes();
ax2.Color = 'none';

yyaxis(ax2,'left')
plot(ax2,Time,N,'LineWidth',1.5)
ylim([1,100])
xlabel('Markov time')
%     ylabel('Number of communities','Color','b')
ylabel('Number of communities')
%     set(ax2(1),'YScale','log','YColor','b')
set(ax2(1),'YScale','log')
set(ax2(1),'XScale','log')

yyaxis(ax2,'right')
plot(ax2,Time,VI,'Color',[0 0.7 0],'LineWidth',1.5)
ylim([0 max(VI)*1.1])
ax2.YAxis(2).Color = [0 0.7 0];
ylabel('Variation of information','Color',[0 0.7 0])

subplot(4,1,[3,4],ax1);
subplot(4,1,[3,4],ax2);

end