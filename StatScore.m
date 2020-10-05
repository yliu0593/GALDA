function [Probability] = StatScore(LDA_xy, Mean_input, std_input)
%StatScore. Will output a probability of the input belonging to a
%particular class. 
%   LDA_xy is a single xy point in LDA Space, MeanLDA is the mean of each
%   class in LDA space, stdLDA is the standard deviation of each class in LDA Space.
%   Written by Youlin Liu & Casey Smith
%   Math by Garth Simpson
%   1/21/2019

%% convert to z-statistic
% z = (x - mu)/std_dev

nClass = length(Mean_input);
for i = 1:nClass
    Z(:,i) = (LDA_xy - Mean_input{i}) ./ std_input{i}; %returns the Z-statistic for LDA_xy belonging to class 1, 2, or 3. 
end

%% Calculate the ratio of probability 
R = zeros(nClass); % preallocate memory for r
Probability = zeros(1,nClass); % preallocate memory for P


%Calculate the values for the r matrix. 
% r is the ratio matrix with values below
    % r11 ... ... ... rn1
    % r12 r22 ... ... ...
    % r13 r33 ... ... ...
    % :
    % :
    % r1n ... ... ... rnn
% Where r12 = exp((-1/2)(z1'*z1 - z2'*z2));
for i = 1:nClass
    for j = 1:nClass
        R(i,j) = exp(-0.5 * (Z(:,i)' * Z(:,i) - Z(:,j)' * Z(:,j))); 
    end   
end
   
%Calculate the Probability, P
    %Where P1 = 1 / sum(R(1,:)), P2 = 1 / sum(R(2,:)), etc. 

for i = 1:nClass
    Probability(i) = 1 / sum(R(:,i));
end


