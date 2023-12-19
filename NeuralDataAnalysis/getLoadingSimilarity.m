function loadingSimilarityVector = getLoadingSimilarity(loadingMatrix)

% based on Umakantha et al., 2021
% inputs:
% loadingMatrix is a matrix of size n x D where n is the number of neurons
% and D is the number of dimensions. 

% outputs:
% a vector of size D where each element is the loading similarity of that
% dimension

% main function
n = size(loadingMatrix,1);
loadingSimilarityVector = 1-(var(loadingMatrix,[],1)./(1/n));