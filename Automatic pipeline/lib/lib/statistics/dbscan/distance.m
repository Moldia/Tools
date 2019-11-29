function d = distance(X,Y, dosqrt)
% function d = distance(X,Y, dosqrt);
% dosqrt: default 0: No square root to reduce computation time (only for comparation perpose)
% Euclidean distance matrix between row vectors in X and Y

if nargin < 3, dosqrt = 0; end

U=~isnan(Y); Y(~U)=0;
V=~isnan(X); X(~V)=0;
if ~dosqrt d=abs(X.^2*U'+V*Y'.^2-2*X*Y');
else
    d=sqrt(X.^2*U'+V*Y'.^2-2*X*Y');
end