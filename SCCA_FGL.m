function [u, v, obj] = SCCA_FGL(X, Y, opts)
% --------------------------------------------------------------------
% Input:
%       - X, e.g., geno matrix
%       - Y, e.g., pheno matrix
%       - opts, parameters: alpha1, alpha2, lambda1, lambda2
% Output:
%       - u, weight of X
%       - v, weight of Y.
%------------------------------------------
% Author: Lei Du, dulei@nwpu.edu.cn
% Date created: 06/05/2017
% Updated: 02/26/2018
%% Copyright (C) 2016-2017 Lei Du
% -----------------------------------------

%% Problem
%
%  max  u' X 'Y v - 1/2*alpha1*||Xu||^2 - 1/2*alpha2*||Yv||^2 -
%  lambda1*||u||_FGL - lambda2*||v||_GGL
% --------------------------------------------------------------------

p = size(X,2);
q = size(Y,2);

% Calculate coverance within X and Y
XX = X'*X;
YY = Y'*Y;

% initialize w1 here
u0 = ones(p, 1);
% --------------------
% initialize w2 here
v0 = ones(q, 1);

u = u0;
v = v0;
% % scale u and v
scale = sqrt(u'*XX*u);
u = u./scale;
scale = sqrt(v'*YY*v);
v = v./scale;

% set parameters
alpha1 = opts.alpha1;
alpha2 = opts.alpha2;
lambda1 = opts.lambda1*opts.lambda1;
lambda2 = opts.lambda2*opts.lambda2;

% get the structure
Gu = updateGraph2(p,'FGL');
Gv = updateGraph2(q,'GGL');

% set stopping criteria
max_Iter = 100;
i = 0;
tol = 1e-5;
obj = [];
tv = inf;
tu = inf;

while (i<max_Iter && (tu>tol || tv>tol)) % default 100 times of iteration
    i = i+1;
    
    % update u
    % -------------------------------------
    % update diagnal matrix D1
    D1 = updateD2(u,Gu,'FGL');
    D1 = diag(D1);
    
    % solve u
    u_old = u;
    F1 = alpha1*XX+lambda1*D1;
    Yv = Y*v;
    b1 = X'*Yv;
    u = F1\b1;
    
    % scale u
    scale = sqrt(u'*XX*u);
    u = u./scale;
    
    % update v
    % -------------------------------------
    % update diagnal matrix D2
    D2 = updateD2(v,Gv,'GGL');
    D2 = diag(D2);
    
    % solve v
    v_old = v;
    F2 = alpha2*YY+lambda2*D2;
    Xu = X*u;
    b2 = Y'*Xu;
    v = F2\b2;
    
    % scale v
    scale = sqrt(v'*YY*v);
    v = v./scale;
    
    % stopping condition
    % -------------------------
    if i > 1
        tu = max(abs(u-u_old));
        tv = max(abs(v-v_old));
    else
        tu = tol*10;
        tv = tol*10;
    end
    obj(end+1) = -u'*X'*Y*v+lambda1*u'*D1*u+lambda2*v'*D2*v+0.5*alpha1*(u'*XX*u-1)+0.5*alpha2*(v'*YY*v-1);
end
end