%Script for performing an adversarial attack. 

%Written by Casey Smith - 12/20/2018
%Updated - 1/23/2019
%Updated - 10/22/2019

%This script is written to assist the LDA Adversarial Spectroscopy project.
%Made to find minimum of the cost function. Original script
%'GeneticAlgorithm_C.m' operated using genetic algorithm to find the
%minimum. However, later work has determined that using the analytical
%solution is much quicker and relatively simple to calculate.

%This script requires the matlab functions "ASolution.m"

clearvars
save_data = 0;

%% User inputs - Import data from LDA program 
numClasses = 3;

% Import Raman Spectra for each class
data{1} = csvread('ProcesseSpec_Class1.csv');
data{2} = csvread('ProcesseSpec_Class2.csv');
data{3} = csvread('ProcesseSpec_Class3.csv');

% Define target class & spectra
initial_class = 3;
target_class = 2;
%spectra_number = round(random('Uniform',1, 84)); %select which spectra you want to attack. 
spectrum_number = 20;

% Import LDA Coordinates
LDAxes = xlsread('LDA_Axes_12262018.xlsx');
LDAxes = LDAxes(:,1:2); %Remove the irrelevant vectors
 LDA_Axes =   LDAxes ;
%%  Perform the Attack

initSpectrum = data{initial_class}(:,spectrum_number);

   %Statistics of spectra 
   for i = 1:numClasses 
       meanSpec{i} = mean(data{i},2);
       varSpec{i} = var(data{i},0,2);
   end

%perturbation = -ASolution(LDAxes, initSpectrum, meanSpec{target_class}); %Don't worry about the negative out front. There was a missed negative sign in the math, it was easier to change here than elsewhere.
perturbation = -ASolution(LDAxes, initSpectrum, data{target_class}(:,20));
perturbedSpectrum = initSpectrum + perturbation; 
   
%%
mask = perturbation;
BaseSpectra = initSpectrum;

% Load all base data
data{1} = csvread('ProcesseSpec_Class1.csv');
data{2} = csvread('ProcesseSpec_Class2.csv');
data{3} = csvread('ProcesseSpec_Class3.csv');

%Statistics of spectra 
for i = 1:3 
   meanSpec{i} = mean(data{i},2);
   stdSpec{i} = std(data{i},0,2); %change to variance? (var)
end
%LDA Data
for i = 1:3
    LDA{i} = transpose(LDA_Axes'*data{i});
    LDA_meanSpec{i} = mean(LDA{i});
    LDA_stdSpec{i} = std(LDA{i});
end

%% Make the Delta

delta = 0:0.01:1;
for i = 1:length(delta)
    delta_mask(:,i) = mask .* delta(i);
    applied_mask(:,i) = delta_mask(:,i) + BaseSpectra;
    scores_mask(:,i) = StatScore(applied_mask(:,i),meanSpec,stdSpec)';
end

delta_LDA = transpose(LDA_Axes'*applied_mask);
for i= 1:length(delta)
    scores_LDA(:,i) = StatScore(delta_LDA(i,:), LDA_meanSpec, LDA_stdSpec)';
end

%% Calculate the magnitude of the perturbation 

perturbation = abs(mask * delta(42));
magnitude = perturbation ./ BaseSpectra; %this will be a percentage

wavenumber = 850:0.47:1479.53;

figure
hold on
title('LDA Space','FontSize',20)
xlabel('LDA 1','FontSize',16)
ylabel('LDA 2','FontSize',16)
plot(LDA{1}(:,1),LDA{1}(:,2),'r.')
plot(LDA{2}(:,1),LDA{2}(:,2),'b.')
plot(LDA{3}(:,1),LDA{3}(:,2),'k.')
plot(delta_LDA(:,1),delta_LDA(:,2),'xg','MarkerSize',8)
plot(delta_LDA(67:87,1),delta_LDA(67:87,2),'xm','MarkerSize',10) %select purple x's, corresponding to the window in plot 3
plot(LDA_meanSpec{1}(:,1),LDA_meanSpec{1}(:,2),'or')
plot(LDA_meanSpec{2}(:,1),LDA_meanSpec{2}(:,2),'ob')
plot(LDA_meanSpec{3}(:,1),LDA_meanSpec{3}(:,2),'ok')
legend('Class 1', 'Class 2', 'Class 3', 'Path of attack', 'Shown in B)','FontSize',14,'Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off


figure
hold on
title('Probability','FontSize',20)
plot(delta, scores_LDA(1,:),'r')
plot(delta, scores_LDA(2,:),'b')
plot(delta, scores_LDA(3,:),'k')
legend('Class 1', 'Class 2', 'Class 3','FontSize',14)
xlabel('|\Delta|','FontSize',16)
ylabel('Probability','FontSize',16)
ax2 = gca;
ax2.FontSize = 16;
hold off

figure
offset1 = 0;
offset2 = 0.002;
offset3 = 0.004;
hold on
title('Initial Spectrum vs. Attacked Spectrum','FontSize',20)
plot(wavenumber, (meanSpec{2}+offset3) / max(meanSpec{2}), 'b')
plot(wavenumber, (applied_mask(:,87)+offset2) / max(meanSpec{2}), 'b') %select appropriate mask
plot(wavenumber, (BaseSpectra+offset1) / max(meanSpec{2}), 'k')
ylabel('Intensity / A.U.','FontSize',16)
xlabel('Wavenumber / cm^{-1}','FontSize',16)
legend('Mean Target Spectrum','Perturbed Spectrum','Initial Spectrum')
ax3 = gca;
ax3.FontSize = 16;
xlim([845, 1500])
ylim([0, 1.8])
hold off



figure
subplot(2,1,1)
plot(wavenumber, (meanSpec{2}) / max(meanSpec{2}), 'k')
title('Mean Target Spectrum')
ylabel('Intensity / A.U.')
%xlabel('Wavenumber / cm^{-1}')
subplot(2,1,2)
plot(wavenumber, (mask * delta(87)) / max(meanSpec{2}), 'm')
title('Applied Perturbation','FontSize',20)
ylabel('Intensity / A.U.','FontSize',16)
xlabel('Wavenumber / cm^{-1}','FontSize',16)
ax4 = gca;
ax4.FontSize = 16;
xlim([845 , 1500])


