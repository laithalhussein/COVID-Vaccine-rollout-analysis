function [vm, va] = extract_vaccine_data_0807_2021a(path_)
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
    if ~ismember(k1, vm_idx) %vaccine administered data
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
            vm.(cmnft).dose2 = 0;
        end
        
        %calculate the total delivered
        vm.(cmnft).dose = vm.(cmnft).dose1 + vm.(cmnft).dose2;
    end    
end

%transform dates into something more tractable
%'Janssen', 'Moderna', 'Pfizer'
va.date_ID = transform_dates(va.date);
vm.Janssen.date_ID = transform_dates(vm.Janssen.date);
vm.Moderna.date_ID = transform_dates(vm.Moderna.date);
vm.Pfizer.date_ID = transform_dates(vm.Pfizer.date);

%% In the vaccine administration data, remove the "Report" data and its corresponding field
report_idx = cellfun(@(x) strcmp(x,'Report'), va.date_type); %14518 instances
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
    disp('Unequal number of territories detected manufacturer data');
    num_vm_terr = NaN;
else
    num_vm_terr = length(vmJ_terr_IDs);
end

%load US state names and abbreviations
load('US_ID_names.mat'); %loads US_full_ID and US_abb_ID
%num_final_US_ID = length(US_full_ID);

%remove the "LTC" territory in the va data (Dont know what that is, and
%not sure what the corollary is in the vm data)
%NOTE: I am assuming that "US" in the va data corresponds to "Federal
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
        keyboard
    end
end

%remove extraneous territories from the vm data
va_terr_IDs = unique(va.terr);
extraneous_terr = vmJ_terr_IDs(find(ismember(vmJ_terr_IDs, va_terr_IDs)==0));

%The extraneous territories are Chicago, NYC, and Philly...it is possible
%that were delievered separately from their respective states. As of now, I
%dunno what to do wit this....

keyboard;

return


function x = transform_dates(dates)
%% transform dates into a numerical value that makes it easy to sort
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