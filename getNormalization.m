function Y = getNormalization(X)
% --------------------------------------------------------------------
% Normalizating data set
% --------------------------------------------------------------------
% Input:
%       - X, input data
% Output:
%       - Y, output data
%------------------------------------------
% Author: Lei Du, leidu@iu.edu
% Date created: DEC-19-2014
% Updated: Jan-16-2015
% @Indiana University School of Medicine.
% -----------------------------------------

[~,n] = size(X);
Y = X;

for j = 1 : n
    Xv = X(:,j);
    Xvn = (Xv-mean(Xv))/std(Xv);
    Y(:,j) = Xvn;
end