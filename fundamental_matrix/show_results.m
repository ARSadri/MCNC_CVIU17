clc
close all
clear

load success_mcnc_once_samples_funda.mat
results_mcnc_1 = myresults;
load success_mcnc_double_samples_funda.mat
results_mcnc_2 = myresults;
load success_mcnc_triple_samples_funda.mat
results_mcnc_3 = myresults;
load success_msnc_once_samples_funda.mat
results_msnc_1 = myresults;
load success_msnc_double_samples_funda.mat
results_msnc_2 = myresults;
load success_msnc_triple_samples_funda.mat
results_msnc_3 = myresults;

figure, hold on
plot(sort(results_mcnc_1), 'r-*')
plot(sort(results_mcnc_2), 'b-*')
plot(sort(results_mcnc_3), 'g-*')
xlim([0 100]); ylim([0 20])
grid on
figure, hold on
plot(sort(results_msnc_1), 'm-*')
plot(sort(results_msnc_2), 'c-*')
plot(sort(results_msnc_3), 'y-*')
xlim([0 100]); ylim([0 20])
grid on
figure, hold on
plot(1, mean(results_mcnc_1, 'all'), 'r*')
plot(2, mean(results_mcnc_2, 'all'), 'b*')
plot(3, mean(results_mcnc_3, 'all'), 'g*')
plot(1, mean(results_msnc_1, 'all'), 'mx')
plot(2, mean(results_msnc_2, 'all'), 'cx')
plot(3, mean(results_msnc_3, 'all'), 'yx')

xlim([0 4]); ylim([0 5])
grid on