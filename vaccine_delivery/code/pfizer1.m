close all;
home;

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultAxesFontSize', 10)
set(0,'defaultTextFontSize', 10)

openfig('pfizer1.fig');
ax = gca;
h = findobj(gca,'Type','line');

x = h.XData;
y = h.YData;
close(gcf);

%convert from categorical array to cell array
xtmp1 = cellstr(x);

%remove the back slashes
for k=1:length(xtmp1), xtmp1{k}([3,6]) = []; end

%cicularly shift the data so that we have a format of year-month-day
xtmp2 = cellfun(@(x) circshift(x,4),xtmp1,'UniformOutput',false);

%resort the data
[xs, sidx] = sort(xtmp2);
ys = y(sidx);

%calculate the instantaneous total of vaccines available for administration
ytotal = cumsum(ys);

%plot it
figure; hold on;
col1 = [0,0,255]/255;
col2 = [255,0,0]/255;
plot([1:length(xs)], ys, 'Linestyle', '-', 'Linewidth', 1.5, 'Marker', '.', 'markersize', 12, 'color', col1);
plot([1:length(xs)], ytotal, 'Linestyle', '-', 'Linewidth', 1.5, 'Marker', '.', 'markersize', 12, 'color', col2);
ylabel('Vaccines delivered');
xlabel(['Weeks aligned to ', cellstr(x(sidx(1)))]);
grid on;

leg1 = legend({'Per week', 'Cumulative distribution'});
set(leg1,'Location', 'Best');


