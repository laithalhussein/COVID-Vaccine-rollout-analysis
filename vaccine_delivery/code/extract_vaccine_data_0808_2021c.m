function [dat, info] = extract_vaccine_data_0808_2021c(path_)
close all;
home;
cd(path_);

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
        va.mavg_scomp = cell2mat(raw_datafile([2:end],15)); %series complete average
        
        va.dose2 = va.dose - va.dose1;
        va.mavg_dose2 = va.mavg_dose - va.mavg_dose1;
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

%create a linear index/key
all_dates = [va.date_num; vm.Janssen.date_num; vm.Moderna.date_num; vm.Pfizer.date_num];
date_seq = sort(unique(all_dates));
%date_key = [1:numel(date_seq)];

%pre-allocate arrays
va.date_ID = nan(size(va.date_num));
for k=1:num_vm
    cmnft = vm_ID{k};
    vm.(cmnft).date_ID = nan(size(vm.(cmnft).date_num));
end

%insert date IDs
for k=1:size(va.date_num,1) %va data
    va.date_ID(k) = find(ismember(date_seq, va.date_num(k)));
end

for k1=1:num_vm %vm data
    cmnft = vm_ID{k1};
    for k2=1:size(vm.(cmnft).date_num,1)
        vm.(cmnft).date_ID(k2) = find(ismember(date_seq, vm.(cmnft).date_num(k2)));
    end
end

%% In the vaccine administration data, remove the "Report" data and its corresponding field
report_idx = cellfun(@(x) strcmp(x,'Report'), va.date_type); %this corresponds to half of all entries
va = rmfield(va, 'date_type');
va_flds = fields(va);
for k=1:length(va_flds), va.(va_flds{k})(report_idx) = []; end

%% Match territories saved in administered and delivered data
%In the vaccine administration data the states are abbreviated. We convert these to the spelled out form
%Note that we also create a territory ID to make things easier to index if need be

%first make sure that we have the same number of territories
vmJ_terr_IDs = unique(vm.Janssen.terr);
vmM_terr_IDs = unique(vm.Moderna.terr);
vmP_terr_IDs = unique(vm.Pfizer.terr);

if length(vmJ_terr_IDs)~=length(vmM_terr_IDs) && length(vmM_terr_IDs)~=length(vmP_terr_IDs)
    disp('Unequal number of territories detected in manufacturer data');
end

%load US state names and abbreviations
load('US_ID_names.mat'); %loads US_full_ID and US_abb_ID
%num_final_US_ID = length(US_full_ID);

%remove the "LTC" territory in the va data (Dont know what that is, and
%not sure what the corollary is in the vm data)
%NOTE: I am also assuming that "US" in the va data corresponds to "Federal
%Entities" in the vm data (check on this!)
LTC_idx = find(cellfun(@(x) strcmp(x, 'LTC'), va.terr));
for k=1:length(va_flds), va.(va_flds{k})(LTC_idx) = []; end

%convert abbreviations
for k=1:size(va.terr,1)
    c_terr = va.terr{k};
    terr_idx = find(cellfun(@(x) strcmp(x, c_terr), US_abb_ID));
    try
        va.terr{k} = US_full_ID{terr_idx};
    catch
        keyboard;
    end
end

%remove extraneous territories from the vm data
va_terr_IDs = unique(va.terr);
extraneous_terr = vmJ_terr_IDs(find(ismember(vmJ_terr_IDs, va_terr_IDs)==0));
disp(['We have ', num2str(length(extraneous_terr)), ' extraneous territories']);

%The extraneous territories are Chicago, NYC, and Philly...it is possible
%that were delievered separately from their respective states. As of now, I
%dunno what to do wit this, but I will add the data to the respective state
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

%% temporally align and window vaccine administered data
%NOTE: We have to assign a date to the data that are windowed.
%To do this, we assign the date that corresponds to the first sample in each window
%If there are discrepnacies in dates between datasets, then we may need to
%shift the onset of this window

%% fill in missing data with zeros (and include corresponding dates)


dat.va = va;
dat.vm = vm;
info = [];

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