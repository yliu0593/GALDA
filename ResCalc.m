function [SquaredRes] = ResCalc(Evecs,Data)
%Inputs: Evecs is an nxm matrix of n evec of horizontal length m
%Data is a class variable of NumSheets classes of input data with each
%column in each class corresponding to a spectrum. 

%This function returns the eigenvalues when the training data are projected
%onto the input eigenvectors. When testing data are used instead, the
%function returns the equivalent projectedd values, calculated analytically
%from projection of the data onto the eigenvectors and using the analytical
%expressions for the Fischer linear discriminant. 

%Determine the number of classes (SheetNum) within the Data
Classes = size(Data);
SheetNum = Classes(1,2);

%Project the data onto the e-vecs
for i=1:SheetNum
    LDA_Projection2{i} = Evecs*Data{i}; %
end

%Calculate the squared resolution for each evec, and return the sum over
%all.
WithinVar=zeros(SheetNum,SheetNum-1); %(number of classes,# of LDA axes)
MeanSpec=zeros(SheetNum,SheetNum-1);
for j = 1:SheetNum-1 %Number of LDA coord 
        for i = 1:SheetNum %number of classes
            WithinVar(i,j) = var(LDA_Projection2{i}(j,:));
            MeanSpec(i,j) = mean(LDA_Projection2{i}(j,:));
         end
        WithinVarLDA(j)=sum(WithinVar(:,j));
        MeanTot(j) = mean(MeanSpec(:,j));
        BtwnVarLDA(j) = (1/SheetNum)*sum((MeanSpec(:,j)-MeanTot(j)).^2);
       %Note - this evaluation implicitly assumes that the number of
       %spectra in each class is identical. If not, the expression needs to
       %be updated to reflect that.
end

EV_calc = BtwnVarLDA./WithinVarLDA; %hell yah! This function recovers the diagonal elements of Vals

SquaredRes = sum(EV_calc);
end