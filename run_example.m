clc;clear;
% -----------------------------------------
% Run this example to see how to use
%------------------------------------------
% Author: Lei Du, dulei@nwpu.edu.cn
% Date created: 02-10-2016
% Update:7-05-2017
% @Northwestern Polytechnical University & Indiana University School of Medicine.
% Contact to: Lei Du or Li Shen (shenli@iu.edu)
% -----------------------------------------

clc; clear;
load example_data.mat

% set parameters, should be tuned before running.
% As an example, we fix them
opts.alpha1 = 1;
opts.alpha2 = 1;
opts.lambda1 = 1;
opts.lambda2 = 0.1;

%% training and testing
[nrow, ~] = size(X);
[test, train] = crossvalind('HoldOut', nrow, 0.7);

X_0 = X(train,:);
Y_0 = Y(train,:);
X_0 = getNormalization(X_0);
Y_0 = getNormalization(Y_0);

XX1 = corr(X');
XX2 = corr(X);
YY1 = corr(Y_0');
YY2 = corr(Y_0);

X_t = X(test,:);
Y_t = Y(test,:);
X_t = getNormalization(X_t);
Y_t = getNormalization(Y_t);

tic;
[u1, v1, obj] = SCCA_FGL(X_0, Y_0, opts);
tt = toc;
corr_XY = corr(X_t*u1,Y_t*v1);

%% results shown
subplot(321)
stem(u0);
title('Ground truth: u');
subplot(322)
stem(v0);
title('Ground truth: v');
subplot(323)
stem(u1);
title('Estimated: u');
subplot(324)
stem(v1);
title('Estimated: v');
subplot(3,2,[5,6])
plot(obj,'-','LineWidth',1.5);
title('Objective value');