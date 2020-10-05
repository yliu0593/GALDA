function [LDA_Projection,Eigenvectors,Eigenvals] = LDAfunc_gjs(Data)
%LDA: Input the data (each class in it's own cell) and function will output
%the eigenvectors.
%   Detailed explanation goes here


SheetNum = numel(Data);
N = zeros(SheetNum, 1); % Number of spectra in each class.
for i=1:SheetNum
   [SpecLength, N(i)] = size(Data{i}); % SpecLength is the number of wavelengths recorded in each spectrum.
end

MeanSpec = zeros(SpecLength, SheetNum); % First index - Wavelengths; Second index - SheetNum. It is a bit counter-intuitive, but it makes sense if you write the entire matrix down.
VarSpec = zeros(SheetNum, SpecLength, SpecLength);
for i=1:SheetNum
    MeanSpec(:, i) = mean(Data{i}, 2);
    VarSpec(i, :, :) = cov(Data{i}'); %generates the within-class variance matrix for each i class, using Matlab built-in covariance function.
end

%Determine the total within-class variance matrix
S_w = squeeze(sum(VarSpec, 1)); % squeeze() function is used here to get rid of the first (SheetNum) dimension in VarSpec.

%Determine the total weighted mean spectrum

Sum_spec_tot=sum(MeanSpec.*repmat(N', SpecLength, 1), 2); % Using A.*repmat(B) can achieve multiplication by row or by column.
Num_spec_tot=sum(N);
Mu_tot = Sum_spec_tot / Num_spec_tot; %this algorithm agrees with MathCad

S_B = zeros(SpecLength);
for i = 1:SheetNum
    S_B = S_B + (MeanSpec(:, i) - Mu_tot)*(MeanSpec(:, i) - Mu_tot)'*(N(i));
end
S_B = S_B/Num_spec_tot;

%DETERMINE THE PROJECTIONS THAT MAXIMIZE THE FISHER LINEAR DISCRIMINANT
%[Vecs,Vals] = eigs(S_B/S_w, SheetNum-1); %note - eig() returns eigenvals as the diagonal
%Note - the above expression is incorrect, yielding S_B*(S_w^-1) rather
%than the correct ordering. 
[Vecs,Vals] = eigs(S_w\S_B, SheetNum-1); %note - eig() returns eigenvals as the diagonal
%[Vecs,Vals] = eigs(inv(S_w)*S_B,SheetNum-1);
%NOTE: Newer Matlab version encourages usage of "Matrix Division" instead of multiplication with inv() when using the inverse as part of a linear algebra expression. This is because inv() function in Matlab is undefined for singular matrices
%elements of a matrix in ascending order, with the most positively valued
%eigenvalues given last. Note - inv(X) returns the inverse matrix of X.

%Next transform the original data by projection onto the linear
%discriminant coordinates. We'll constrain the plotting to the first two
%(assuming at least three classes of data).

%Project the raw data onto the selected principal components within
%'transform' and write to an output file

for i=1:SheetNum
    LDA_subset{i}= Vecs'*Data{i};
    LDA_subset_trans{i}=LDA_subset{i}'; %Transpose operation arranges the 
    %data as XY pairs for 2D plotting and XYZ triplets for 3D plotting.
end

LDA_Projection = LDA_subset_trans;
Eigenvectors = real(Vecs');
Eigenvals = real(Vals);
end

