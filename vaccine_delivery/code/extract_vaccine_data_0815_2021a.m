function [dat] = extract_vaccine_data_0815_2021a(path_, LTC_flag, exclude_onset_flag)
close all;
home;
cd(path_);

%define flag which determines whether or not we want to exclude the "LTC"
%territory from the adminstrated dataset
if nargin<2
    LTC_flag = 0;
end

%define flag which determines if we want to add the administered vaccines
%from the very day recorded with its following week (necessary for temporal
%alignment)
if nargin<3
    exclude_onset_flag = 0;
end

dbstop if caught error

%% extract the raw data first

%define vaccine manufacturers relevant for US 
%(needs to case-match the name in the file)
vm_ID = {'Janssen', 'Moderna', 'Pfizer'};
num_vm = length(vm_ID);

%assume ALL data are stored in CSV files
d=dir([path_,'\*.csv']);

%save filenames
count = 0;
filename_list = cell(1,length(d));
for k=1:length(d)
    if d(k).isdir==0
        count=count+1;
        filename_list{count}=d(k).name;
    end
end
num_files = count;

%determine the index of which files correspond to which vaccine manufacturers
vm_idx = nan(1,num_vm);
cc = 1;
for k1=1:num_files    
    for k2=1:num_vm
        if contains(filename_list{k1}, vm_ID{k2})
           vm_idx(cc) = k1;
           cc=cc+1;           
        end
    end
end

%determine index of vaccine administered data (assume there's only one)
va_idx = find(ismember([1:num_files], vm_idx)==0);

%set up some flags to determine if we're dealing with datasets that have
%state data
va_national_flag = 1;
vm_national_flag = 1;

%save the data
for k1=1:num_files
    cfile = filename_list{k1};
    [~,~,raw_datafile] = xlsread(cfile);
    
    %tailor extraction of the data depending on what we're working with
    if k1==va_idx %vaccine administered data
        va.date = raw_datafile([2:end],1);
        va.terr = raw_datafile([2:end],3);
        va.dose = cell2mat(raw_datafile([2:end],4)); %daily administration
        va.mavg_dose = cell2mat(raw_datafile([2:end],6)); %moving average administration
        va.dose1 = cell2mat(raw_datafile([2:end],7)); %daily scale
        va.mavg_dose1 = cell2mat(raw_datafile([2:end],9));
        va.date_type = raw_datafile([2:end],10); %takes either "Admin" or "Report" (only trust Admin)
        va.scomp = cell2mat(raw_datafile([2:end],13)); %series complete
        va.mavg_scomp = cell2mat(raw_datafile([2:end],15)); %series complete moving average
        
        va.dose2 = va.dose - va.dose1;
        va.mavg_dose2 = va.mavg_dose - va.mavg_dose1;
        
        %determine if we have local territory data
        if length(unique(va.terr))>=50
            va_national_flag = 0;
        end
        
    else %vaccine manufacturer data
        cmnft = vm_ID{find(k1==vm_idx)};
        vm.(cmnft).terr = raw_datafile([2:end],1);
        vm.(cmnft).date = raw_datafile([2:end],2);
        vm.(cmnft).dose1 = cell2mat(raw_datafile([2:end],3));
        
        if size(raw_datafile,2)>3
            vm.(cmnft).dose2 = cell2mat(raw_datafile([2:end],4));
        else
            vm.(cmnft).dose2 = vm.(cmnft).dose1*0;
        end
        
        %calculate the total delivered
        vm.(cmnft).dose = vm.(cmnft).dose1 + vm.(cmnft).dose2;
        
        %determine if we have local territory data for ALL manufacturer data
        if length(unique(vm.(cmnft).terr))>=50
            vm_national_flag = 0;
        end
    end    
end

%% dates are entered as XX/XX/XXXX or X/XX/XXXX or XX/X/XXXX - I standardize this below
va.date = standardize_dates(va.date);
vm.Janssen.date = standardize_dates(vm.Janssen.date);
vm.Moderna.date = standardize_dates(vm.Moderna.date);
vm.Pfizer.date = standardize_dates(vm.Pfizer.date);

%% create a date ID that will permit linear indexing wrt time across ALL datasets
%'Janssen', 'Moderna', 'Pfizer'
va.date_num = transform_dates(va.date);
vm.Janssen.date_num = transform_dates(vm.Janssen.date);
vm.Moderna.date_num = transform_dates(vm.Moderna.date);
vm.Pfizer.date_num = transform_dates(vm.Pfizer.date);

%% In the vaccine administration data, remove the "Report" data and its corresponding field
report_idx = cellfun(@(x) strcmp(x,'Report'), va.date_type); %this corresponds to half of all entries
va = rmfield(va, 'date_type');
va_flds = fields(va);
for k=1:length(va_flds), va.(va_flds{k})(report_idx) = []; end

%% Match territories saved in administered and delivered data
%In the vaccine administration data the states are abbreviated. We convert these to the spelled out form
%Note that we also create a territory ID to make things easier to index if need be

%first make sure that we have the same number of territories for each
%manufacturer
if vm_national_flag==0
    vmJ_terr_IDs = unique(vm.Janssen.terr);
    vmM_terr_IDs = unique(vm.Moderna.terr);
    vmP_terr_IDs = unique(vm.Pfizer.terr);
    
    if length(vmJ_terr_IDs)~=length(vmM_terr_IDs) && length(vmM_terr_IDs)~=length(vmP_terr_IDs)
        disp('Unequal number of territories detected in manufacturer data');
    end
end

%optionally remove the "LTC" territory in the va data
if LTC_flag
    LTC_idx = find(cellfun(@(x) strcmp(x, 'LTC'), va.terr));
    for k=1:length(va_flds), va.(va_flds{k})(LTC_idx) = []; end
end

%convert abbreviations if we're dealing with state va data
if va_national_flag==0
    %load US state names and abbreviations
    [US_full_ID, US_abb_ID] = extract_suppl_data;
    
    for k=1:size(va.terr,1)
        c_terr = va.terr{k};
        terr_idx = find(cellfun(@(x) strcmp(x, c_terr), US_abb_ID));
        try
            va.terr{k} = US_full_ID{terr_idx};
        catch
            keyboard;
        end
    end
end

%remove extraneous territories from the vm data
if ~va_national_flag && ~vm_national_flag
    va_terr_IDs = unique(va.terr);
    extraneous_terr = vmJ_terr_IDs(find(ismember(vmJ_terr_IDs, va_terr_IDs)==0));
    disp(['We have ', num2str(length(extraneous_terr)), ' extraneous territories']);
end

%In the vm data, we have extraneous territories in Chicago, NYC, and Philly
%it is possible that were delievered separately from their respective states.
%As of now, I dunno what to do with this, but I will add the data to the respective state
if vm_national_flag
    vm_flds = fields(vm.Pfizer);
    for k1=1:num_vm
        cmnft = vm_ID{k1};
        
        IL_idx = find(cellfun(@(x) strcmp(x,'Illinois'), vm.(cmnft).terr));
        Chicago_idx = find(cellfun(@(x) strcmp(x,'Chicago'), vm.(cmnft).terr));
        
        NY_idx = find(cellfun(@(x) strcmp(x,'New York'), vm.(cmnft).terr));
        NYC_idx = find(cellfun(@(x) strcmp(x,'New York City'), vm.(cmnft).terr));
        
        PA_idx = find(cellfun(@(x) strcmp(x,'Pennsylvania'), vm.(cmnft).terr));
        Philadelphia_idx = find(cellfun(@(x) strcmp(x,'Philadelphia'), vm.(cmnft).terr));
        
        %combine city data with their respective states
        vm.(cmnft).dose1(IL_idx) = vm.(cmnft).dose1(IL_idx) + vm.(cmnft).dose1(Chicago_idx);
        vm.(cmnft).dose1(NY_idx) = vm.(cmnft).dose1(NY_idx) + vm.(cmnft).dose1(NYC_idx);
        vm.(cmnft).dose1(PA_idx) = vm.(cmnft).dose1(PA_idx) + vm.(cmnft).dose1(Philadelphia_idx);
        
        vm.(cmnft).dose2(IL_idx) = vm.(cmnft).dose2(IL_idx) + vm.(cmnft).dose2(Chicago_idx);
        vm.(cmnft).dose2(NY_idx) = vm.(cmnft).dose2(NY_idx) + vm.(cmnft).dose2(NYC_idx);
        vm.(cmnft).dose2(PA_idx) = vm.(cmnft).dose2(PA_idx) + vm.(cmnft).dose2(Philadelphia_idx);
        
        %remove the cities
        for k2=1:length(vm_flds)
            cvm_fld = vm_flds{k2};
            vm.(cmnft).(cvm_fld)([Chicago_idx, NYC_idx, Philadelphia_idx]) = [];
        end
    end
end

if ~va_national_flag && ~vm_national_flag
    tflag = 0;
    for k=1:length(vm_ID)
        cmnft = vm_ID{k};
        if length(unique(va.terr)) ~= length(unique(vm.(cmnft).terr)), tflag = 1; end
    end
    if tflag==1, warning('We have territory discontinuity between datasets!'); end
end

%% temporally align and window vaccine administered data
%NOTE: We have to assign a date to the data that are windowed.
%To do this, we assign the date that corresponds to the first sample in each window
%If there are discrepnacies in dates between datasets, then we may need to
%shift the onset of this window

wsize = 7;

%we will create a 3D matrix for the windowed va data with dimensions:
%week x territory x dose
total_weekly_dates = ceil((numel(unique([va.date_num]))-1)/wsize); %subtract by one because we dont count the first day in the va data
total_va_territories = numel(unique(va.terr));

wva_data = zeros(total_weekly_dates, total_va_territories, 2);
win = [2:wsize:numel(unique(va.date_num))]; %indexing from 2 is for ignoring the first day
    
%loop through each territory
for k1=1:numel(unique(va.terr))
    
    %find all indices that pertain to the current territory
    cterr_idx = find(cellfun(@(x) strcmp(x, va.terr{k1}), va.terr));
    
    %acquire the corresponding dates and sort them
    [cdates_num, cdates_num_idx] = sort(va.date_num(cterr_idx));
    
    %IMPORTANT: There is a discrepnacy between the very first day the vaccine was
    %administered (12/13/2020) and the firt day it was delivered (12/14/2020) 
    %Exclusion of this data may not be the best strategy, so I optionally 
    %add the administered from this day onto the adminitered for the following week
    
    %we window the data in the temporal order in which the vaccines were administered
    wdose1 = window_data_V2(va.dose1(cterr_idx(cdates_num_idx(2:end))),wsize,@nansum);
    wdose2 = window_data_V2(va.dose2(cterr_idx(cdates_num_idx(2:end))),wsize,@nansum);
    
    if exclude_onset_flag==0 %add first day onto the first week
        %Note: this addition is naive if we wanted to combine the data in
        %any other way, like taking the mean for example
        wdose1(1) = wdose1(1)+va.dose1(cterr_idx(cdates_num_idx(1)));
        wdose2(2) = wdose2(1)+va.dose2(cterr_idx(cdates_num_idx(1)));
    end
    
    %save it
    wva_data(:,k1,1) = wdose1;
    wva_data(:,k1,2) = wdose2;
    
    %acquire and save the corresponding dates for these weeks, which are
    %identical for all states
    if k1==1
        wdates_num = cdates_num(win);
        wdate_tmp = va.date(cdates_num_idx);
        wdate = wdate_tmp(win);
    end
end

%save this windowed data into a new structure
V.a = wva_data;
V.date_num = wdates_num;
V.date = wdate;

%% ensure that vaccine manufacturer data are analogously aligned

%NOTE: I assume that we always have more up-to-date administered data
%compared to manufacturer data (i.e., the span is larger) because it is
%sampled daily, not weekly

%make new manufacturer data organized by the following:
%date x territory x manufacturer (J/M/P) x dose allocation
total_vm_territories = numel(unique(vm.(cmnft).terr));

%Note that we still temporally align the datasets even if we dont match the
%terriory data
if ~vm_national_flag && ~va_national_flag
    wvm_data = zeros(total_weekly_dates, total_va_territories, 3, 2);
else
    wvm_data = zeros(total_weekly_dates, total_vm_territories, 3, 2);
end

%If we have state data for both va and vm, then we order the state indices
%of the vm data in the same order of territories as that of the va data
if ~vm_national_flag && ~va_national_flag
    for km = 1:length(vm_ID) %loop through each manufacturer
        cmnft = vm_ID{km};
        for k1=1:numel(unique(va.terr)) %index the same way we did through administration data
            
            cterr_idx = find(cellfun(@(x) strcmp(x, va.terr{k1}), vm.(cmnft).terr)); %note we already made sure these match in content
            [cdates_num, original_dates_idx] = sort(vm.(cmnft).date_num(cterr_idx)); %find when in time this territory acquired the vaccine
            
            %find where these dates equate to the corresponding date indices of
            %the administered data
            new_dates_idx = find(ismember(V.date_num, cdates_num));
            
            %insert the data in their corresponding indices
            wvm_data([new_dates_idx], k1, km, 1) = vm.(cmnft).dose1(cterr_idx(original_dates_idx));
            wvm_data([new_dates_idx], k1, km, 2) = vm.(cmnft).dose2(cterr_idx(original_dates_idx));
        end
    end
else
    for km = 1:length(vm_ID)
        cmnft = vm_ID{km};
        for k1=1:total_vm_territories
            cterr_idx = find(cellfun(@(x) strcmp(x, vm.(cmnft).terr{k1}), vm.(cmnft).terr));
            [cdates_num, original_dates_idx] = sort(vm.(cmnft).date_num(cterr_idx));
            new_dates_idx = find(ismember(V.date_num, cdates_num));
            wvm_data([new_dates_idx], k1, km, 1) = vm.(cmnft).dose1(cterr_idx(original_dates_idx));
            wvm_data([new_dates_idx], k1, km, 2) = vm.(cmnft).dose2(cterr_idx(original_dates_idx));
        end
    end
end

V.m = wvm_data;
dat.va = va;
dat.vm = vm;
dat.V = V;
return

function x = transform_dates(dates)
%transform dates into a numerical value that makes it easy to sort
%dates are in month/day/year format

N = size(dates,1);

%remove the back slashes in vaccine administered data
date_tmp = dates;
date_ID_tmp = cell(N,1);
for k=1:N
    date_k_tmp = date_tmp{k};
    date_k_tmp(date_k_tmp=='/') = [];
    
    %cicularly shift the date ID so that we have a format of year-month-day
    %(important so that we can sort the dates correctly)
    date_ID_tmp{k} = str2double(circshift(date_k_tmp,4));
end

x = cell2mat(date_ID_tmp);
return

function x = standardize_dates(dates)
%make all dates so that they are in XX/XX/XXXX format

N = size(dates,1);
new_date_format = cell(N,1);

for k=1:N
    date_k_tmp = dates{k};
    
    %determine location of the '/' symbol
    marker_idx = find(date_k_tmp=='/');
    string_A = date_k_tmp(1:marker_idx(1)-1);
    string_B = date_k_tmp(marker_idx(1)+1:marker_idx(2)-1);
    string_C = date_k_tmp(marker_idx(2)+1:end);
    
    string_format = [marker_idx(1), diff(marker_idx)]-1;
    if string_format(1)==1
        new_string_A = ['0',string_A];
    else
        new_string_A = string_A;
    end
    if string_format(2)==1
        new_string_B = ['0',string_B];
    else
        new_string_B = string_B;
    end
    
    new_date_format{k} = [new_string_A, '/', new_string_B, '/', string_C];
end

x = new_date_format;
return