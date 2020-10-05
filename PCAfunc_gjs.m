function [Projection,Eigenvectors] = PCAfunc_gjs(Data, OutputDim)
%PCA: Input the data (each class in it's own cell) and function will output
%the eigenvectors.
%   Detailed explanation goes here
% OutputDim is the dimension of the desired PCA vectors returned by the function.


SheetNum = numel(Data);
N = zeros(SheetNum, 1); % Number of spectra in each class.
for i=1:SheetNum
   [SpecLength, N(i)] = size(Data{i}); % SpecLength is the number of wavelengths recorded in each spectrum.
end

for i=1:SheetNum
    if i == 1
        DataPooled = Data{i};
    else
        DataPooled = horzcat(DataPooled,Data{i});
    end
end

%DETERMINE THE PRINCIPAL COMPONENTS
%Find the mean vector. The number of columns is given by the second entry
%returned from the function 'size'.

DataSize=size(DataPooled);
for i=1:DataSize(2) %sum over the number of columns in the pooled dataset
    if i == 1
      SummedVec=DataPooled(:,i); %initialize by setting equal to the first column
    else
       SummedVec=SummedVec+DataPooled(:,i);
    end
end

MeanVec = (1/DataSize(2))*SummedVec;
%zero-center the spectra
for i=1:DataSize(2) %Evaluate on a column by column basis
    DataPooledZC(:,i)=DataPooled(:,i)-MeanVec; 
end
%Generate the variance covariance matrix alpha by mult by the transpose
alpha = DataPooledZC*DataPooledZC.';
[Vecs,Vals] = eig(alpha); %note - eig() returns eigenvals as the diagonal 
%elements of a matrix in ascending order, with the most positively valued
%eigenvalues given last. 

%Next transform the original data by projection onto the principal
%component vectors given by the eigenvecs. Note - because the eigenvalues
%are ranked in ascending order, the last eigenvectors are the relevant
%ones, working backwards.
for i=0:OutputDim-1
    if i == 0
        transform(:,1)=Vecs(:,length(Vecs)-i);
    else
        transform = horzcat(transform,Vecs(:,length(Vecs)-i));
    end
end

%Project the raw data onto the selected principal components within
%'transform'

for i=1:SheetNum
    PC_subset{i}= transform'*Data{i};
    PCsubset_trans{i}=PC_subset{i}'; %Transpose operation arranges the 
    %data as XY pairs for 2D plotting and XYZ triplets for 3D plotting.
end

%Projection = PCsubset_trans;
Projection = PC_subset; %note - there are times to use the transpose as well. 
%When coupling with LDA, the straight-up projection must be used.

Eigenvectors = real(transform');
end

