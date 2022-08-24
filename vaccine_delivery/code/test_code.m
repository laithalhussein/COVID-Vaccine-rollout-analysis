clear all;
close all;
home;

my_path = 'C:\Users\laith\Dropbox (HNL)\Laith_files\vaccine_delivery\data';

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultAxesFontSize', 10)
set(0,'defaultTextFontSize', 10)

%acquire data
[dat] = extract_vaccine_data_0814_2021b(my_path);
