function [d, e, f] = fit_circle_nhom(X)
% function [d e f] = fit_circle_nhom(X)
%
% INPUT: 
% X: n-by-2 vector
%    with data
%
%
% OUTPUT: 
% quadric coordinates of the circle
N = size(X,1);
I = ones(N,1); % sloupec jedniček
A = [X I];
b = -sum(X.^2, 2);

params = A \ b;
d = params(1);
e = params(2);
f = params(3);

end
