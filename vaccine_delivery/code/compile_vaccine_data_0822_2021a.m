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
vr = cumsum(vmtb) - cumsum(vatb);
plot(t, vr, 'color', 'b', 'linewidth', 1.5);

%plot a vertical line where vaccines where available for all adults
ta_idx = find(V.date_num==20210419);
ta = t(ta_idx);
vline(ta, 'color', 'k', 'linestyle', '--');

leg_tmp = legend({'Vaccines delivered: V_d(t)', 'Vaccines administered: V_a(t)',...
    'Vaccines remaining: V_R(t) = \SigmaV_d(t) - \SigmaV_a(t)', 'Week vaccines were available to all'});
set(leg_tmp, 'Location', 'Best');

ylabel('Number of vaccines');
xlabel('Weeks since vaccine availability');
title('Time series of vaccine activity');
grid on;

%% calculate metrics for demand

%Conservative metric: D(t) = V_a(t) + alpha*V_a(t+1) (note that demand will
%always be greater than V_a in this case)
va_offset = [vatb(2:end); 0];
Na = 20;
alpha = linspace(0,1,Na+1); %use a range (0-100%) of alpha values

%first calculate alpha*V_a(t+1)
dtmp = (alpha'*va_offset')';

%add V_a(t) onto each column of dtmp to caluclate demand
d=nan*dtmp;
for k=1:Na+1, d(:,k) = dtmp(:,k)+vatb(k); end

%plot vaccines remaining with demand for 3 casesL d = 0, d = 0.5, and d = 1
figure; hold on;
plot(t, vr, 'color', 'b', 'linewidth', 1.5);
plot(t, vatb, 'color', brown, 'linewidth', 1.5);
plot(t, d(:,1), 'color', yellow, 'linestyle', '--', 'linewidth', 1.5);
plot(t, d(:,floor(Na/2)), 'color', orange, 'linestyle', '--', 'linewidth', 1.5);
plot(t, d(:,end), 'color', pink, 'linestyle', '--', 'linewidth', 1.5);

leg_tmp = legend({'Vaccines remaining: V_R(t)', 'Vaccines administered', 'D(t) for \alpha = 0',...
    'D(t) for \alpha = 0.5', 'D(t) for \alpha = 1'});
set(leg_tmp, 'Location', 'Best');

ylabel('Number of vaccines  / Estimate of demand for vaccines');
xlabel('Weeks since vaccine availability');
title('Vaccines remaining compared to Demand (D(t) = V_R(t) + \alphaV_R(t+1))');

%now calculate fraction of vaccines remaining that are not due to demand
vr2 = nan*d;
for k=1:Na+1, vr2(:,k) = d(:,k)-vatb(k); end

TI = nan*vr2; %throughput index
for k=1:Na+1, TI(:,k) = vr2(:,k)/vr(k); end

%visualize the result
h1 = figure;
imagesc(TI');

yt = get(gca, 'YTick');
set(gca, 'YTickLabel', alpha);
set(gca, 'YTick', yt, 'YTickLabel', alpha(yt), 'YTickLabelRotation', 0);

ylabel('Level of demand (\alpha)');
xlabel('Weeks since vaccine availability');
title('Fraction of vaccines remaining not due to demand');

hc = colorbar;
ylabel(hc, 'Throughput Index (TI)');


%Note between weeks 4-18, slope of V_R is negative, so demand here is
%unquestionable, and we should focus analysis here
%%


%keyboard;