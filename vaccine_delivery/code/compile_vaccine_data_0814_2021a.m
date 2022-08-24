clear all;
close all;
home;

my_path = 'C:\Users\laith\Dropbox (HNL)\Laith_files\vaccine_delivery\data';

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultAxesFontSize', 10)
set(0,'defaultTextFontSize', 10)

purple = [0.5,0,0.5];
grey = [0.5,0.5,0.5];
orange = [255, 165, 0]/255;
dgreen = [34, 139, 34]/255;
brown = [210,105,30]/255;
light_grey = [211,211,211]/255;
cyan = [0,255,255]/255;
pink = [255, 105, 180]/255;
yellow = [255, 179, 0]/255;

%acquire data
%manufacturer data dims: date x territory x manufacturer (J/M/P) x dose allocation
%administration data dims: x territory x dose #
[dat] = extract_vaccine_data_0815_2021a(my_path);
[V, va, vm] = deal(dat.V, dat.va, dat.vm);
%% plot instantaneous vaccines delieverd and administered

%define time vector
t = [1:length(V.date_num)]-1; %in weeks

%calculate total vaccines delieverd first across manufacturers and doses
vmta = squeeze(sum(nansum(V.m,3),4));
vata = squeeze(nansum(V.a,3));

%sum across territories
vmtb = sum(vmta,2);
vatb = sum(vata,2);

figure; hold on;
plot(t, vmtb, 'color', cyan, 'linewidth', 1.5);
plot(t, vatb, 'color', brown, 'linewidth', 1.5);

%calculate and plot total vaccines available based on the amount delievered
%on a given week and the amount available from previous weeks
v_av = cumsum(vmtb) - cumsum(vatb);
plot(t, v_av, 'color', 'b', 'linewidth', 1.5);

%plot a vertical line where vaccines where available for all adults
ta_idx = find(V.date_num==20210419);
ta = t(ta_idx);
vline(ta, 'color', 'k', 'linestyle', '--');

leg_tmp = legend({'Vaccines delivered: V_d(t)', 'Vaccines administered: V_a(t)',...
    'Vaccines remaining: \SigmaV_d(t) - \SigmaV_a(t)', 'Week vaccines were available to all'});
set(leg_tmp, 'Location', 'Best');

ylabel('Number of vaccines');
xlabel('Weeks since vaccine availability');
title('Time series of vaccine activity');
grid on;

%% calculate a metric for throughput that accounts for demand

%NOTES:
%nominally, a difference between va and vm indicates a discrepancy, but how
%would we know if this is due to weak throughput or lack of demand?

%For any given point in time, the "ideal" amount of vaccines to be
%administered should be equal to the amount of vaccines that are available,
%assuming that demand meets or exceeds the number of vaccines available. If
%demand were to be lower than the vaccines that are available, then the
%"ideal" amount of vaccines to administer should be equal to demand
%instead. Thus we need an index that quantifies how far away we are from
%this "ideal"

%%


keyboard;