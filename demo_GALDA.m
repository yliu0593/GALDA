%% Generative Adversarial Spectroscopy

%Written by Casey Smith 
    %1st edition - 5/14/2019
    %Updated - 11/15/2019: Target changed from mean of target class, to
    %individual spectrum of target class. See line 129.
    %Updated - 11/19/2019: Updated LDAfunc to LDAfunc2 which now outputs
        %eigenvalues too
    %Updated on 5/31/2020 by Garth Simpson to divide data into testing
        %and training sets, and to provide statistics on overfitting as a
        %function of attack progression.

        
%Overview:
    %Perform LDA on a data set
    %Generate random seeds (junk data) 
    %Cause junk data to misclassify to various classes
    %Perform LDA again, but now with the junk data as an additional class
    %Iterate (until impatience?)
    
clearvars;
AdditionalSpec = 0; % Number of additional spectra added to the experimental data set to support matrix inversion
TestingFraction = 0.2; % 0.2(67); 0.8(16), 0.9(8), 0.96(3) %Fraction of the genuine data held out for testing
TrainingFraction = 1-TestingFraction;
NumGenSpec =100; %Number of generated adversarial spectra in each iterative loop
NumLoops = 80;
ShuffleData = 1;

%% Import data set
%This script is designed for spectra to be saved in an excel file, with
%each classes saved on different sheets.

%Retrieve info on sheets (names & number)
InputFilename = 'ProcesseSpec_All.xlsx';
[~, InputFileSheets] = xlsfinfo(InputFilename); % Get info about input spread sheet for automatic class counting. 
SheetNum = numel(InputFileSheets); % Automatically determine the number of worksheets in input file.


for i = 1:SheetNum
    data_tot_beforeShuffle{i} = xlsread(InputFilename, InputFileSheets{i}); %Save each sheet in a different cell of the variable 'data'
    datasize_tot(i,:) = size(data_tot_beforeShuffle{i});
    if ShuffleData == 1
        data_tot{i} = data_tot_beforeShuffle{i}(:,randperm(datasize_tot(i,2)));
    else
        data_tot{i} = data_tot_beforeShuffle{i};
   end
    NumTraining = floor(TrainingFraction*datasize_tot(i,2));
    data{i}= data_tot{i}(:,1:NumTraining);
    dataTest{i} = data_tot{i}(:,NumTraining:datasize_tot(i,2));
    datasize(i,:) = size(data{i});
    datasizeTest(i,:) = size(dataTest{i});
end

%% LDA Prep
%We need to generate enough additional spectra to make the system
%solvable. E.g. 84 spectra with 1340 channels would require at
%least a couple hundred additiona spectra to be evaluable.


for i=1:SheetNum
    meanSpec{i} = mean(data{i},2); %determine the mean spectra of each class
    stdSpec{i} = std(data{i},0,2); %determine the standard deviation at each wavelength for each class
    sizeInfo{i} = size(data{i}); %save the matrix dimensions of each class to determine how many spectra need to be generated.
    SpecLength = max(sizeInfo{i});
end

%generate new spectra based off of statistics
for i=1:SheetNum
    SpecOut = zeros(length(data{i}(:,1)),AdditionalSpec);
    %for j=1:max(sizeInfo{i})
    for j=1:AdditionalSpec
        New_Spec = normrnd(meanSpec{i}, stdSpec{i});
        SpecOut(:,j) = New_Spec;
    end
    GeneratedData{i} = horzcat(data{i},SpecOut);  
    clear SpecOut New_Spec j k
end

%% Perform Initial LDA 

[LDA_Projection, Vecs, Vals{1}] = LDAfunc_gjs(GeneratedData); %LDA_Projection is all of the data (including generated data)
Vecs(3,:) = 0;
%Project testing and training data onto the eigenvecs
for i=1:SheetNum
    LDA_Projection2{i} = Vecs*data{i}; %Projection2 is just the genuine training data
    LDA_Projection2{i} = LDA_Projection2{i}';
    LDA_Projection_Testing{i} = Vecs*dataTest{i};
    LDA_Projection_Testing{i} = LDA_Projection_Testing{i}';
end

%Recover the eigenvalues for the original training and testing data
EVSum(1) = ResCalc(Vecs,data);
EVSumTesting(1) =ResCalc(Vecs,dataTest);
Vecs(3,:) = 0;
%perform projection of just original data
clr{1}=[1 0 0];
clr{2}=[0 0 1];
clr{3}=[0 0 0];
clr{4}=[1 0 1];



for i = 1:SheetNum
    SeedLowerBound(i) = min(min(data{i})); %determine the lowest value for the uniform distribution
    SeedUpperBound(i) = max(max(data{i})); %determine the highest value for the uniform distribution
end


%% Generate Random Seeds
loop = 1;
while loop < NumLoops

SavedVecs{loop} = Vecs;

SeedData = random('uniform', min(SeedLowerBound)*.7, max(SeedUpperBound)*.7, [SpecLength,NumGenSpec]);
SeedLDA = transpose(Vecs*SeedData);

%% Generate Decoy data to various classes

targetClass = randi(3,[size(SeedData,2), 1]);
targetSpec = randi(min(min(datasize)), [size(SeedData,2), 1]);

GenStatus = 1;
parfor i = 1:size(SeedData,2)
    %JunkData_Perturbed(:,i) = GA_Func(JunkData(:,i),targetClass(i),Vecs, meanSpec);
    DecoyData_Perturbed(:,i) = SeedData(:,i) - ASolution(Vecs', SeedData(:,i), meanSpec{targetClass(i)});
    %JunkData_Perturbed(:,i) = JunkData(:,i) - ASolution(Vecs', JunkData(:,i), data{targetClass(i)}(:,targetSpec(i)));
    GenStatus = GenStatus + 1
end

if loop == 1
    GeneratedData{4} = DecoyData_Perturbed;
else
    GeneratedData{4} = horzcat(GeneratedData{4},DecoyData_Perturbed);
end

%perform LDA again
[LDA_Projection, Vecs, Vals{loop}] = LDAfunc_gjs(GeneratedData);

%%Project the genuine data onto just the relevant coords of vecs (i.e.,
%remove first column, cooresponding to decoy data direction). 
Vecs2 = Vecs(2:SheetNum,:);
%Next project testing and training data onto the relevant eigenvecs
for i=1:SheetNum
    LDA_Projection2{i} = Vecs2*data{i}; 
    LDA_Projection2{i} = LDA_Projection2{i}';
    LDA_Projection_Testing{i} = Vecs2*dataTest{i};
    LDA_Projection_Testing{i} = LDA_Projection_Testing{i}';
end

%Recover the Evals for the training and testing data post-attack
EVSum(loop+1) = ResCalc(Vecs2,data);
EVSumTesting(loop+1) = ResCalc(Vecs2,dataTest);

loop = loop+1
end



EVSum = abs(EVSum');
EVSumTesting = abs(EVSumTesting');



%%
LDA_Projection2{1} = (Vecs*data{1})';
LDA_Projection2{2} = (Vecs*data{2})';
LDA_Projection2{3} = (Vecs*data{3})'; 
LDA_decoy = (Vecs*GeneratedData{4})';
set(0,'DefaultFigureWindowStyle','normal')
figure1 = figure('units','normalized','outerposition',[0 0 0.8 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
for i=1:SheetNum
   scatter3(LDA_Projection2{i}(:,3),LDA_Projection2{i}(:,2),LDA_Projection2{i}(:,1),14,clr{i},'filled');
end
colormap(axes1,spring)
scatter3(LDA_decoy(:,3),LDA_decoy(:,2),LDA_decoy(:,1),70,1:NumGenSpec*(loop-1),'x');
xlabel('LD3');
ylabel('LD2');
zlabel('LD1');
set(axes1,'FontSize',20);
colorbar
legend('Class 1', 'Class 2', 'Class 3','Decoy Data')
