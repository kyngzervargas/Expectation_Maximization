clear all;
% close all;

load 'apple(hue,ecc).mat'   
load 'banana(hue,ecc).mat'
load 'orange(hue,ecc).mat'
load labels.mat

% YY = readtable('DataYY.xlsx');

%% 
% no_ = ones(1,size(apple,1))*-1;
% % no_ = zeros(1,size(apple,1)); %FOR ACT 14
% % no_ = ones(1,size(banana,1))*-1;
% o_ = ones(1,size(orange,1));
% no_ = no_';
% o_ = o_';
% 
% %LABEL matrix
% d = [o_;no_];
% % d = flipud(d);
% 
% %Data matrix
data = [apple;orange];
% % data = [banana;orange];
% 
% %Variables
% bias = ones(1,size(data,1));
% bias = bias';
% n = 0.5; %learning rate
% 
% X = [bias,data]; %Data
% 
% %WEIGHTS
% a_ = 0; b_ = 1; %range
% w = (b_-a_).*rand(3,1)/100 + a_;
% w = w';
% % weights = -1*2.*rand(3,1);

%%

% A = ['orange'];
% X = 20;
% Output1 = repmat({A},X,1);
% Y = [Output;Output1];

XX = data;
YY = labels;
% YY = {'apple','orange'};
% YY = [
% YY = [['apple';'orange'],[20;20],[50;50]];


%%

Mdl = fitcnb(XX,YY, 'ClassNames',{'apple','orange'});
setosaIndex = strcmp(Mdl.ClassNames,'apple');
estimates = Mdl.DistributionParameters{setosaIndex,1};

%%
CLR = ['b','k'];
SYM = 'xd';
figure
gscatter(XX(:,1),XX(:,2),YY,CLR,SYM);
h = gca;
cxlim = h.XLim;
cylim = h.YLim;
hold on
Params = cell2mat(Mdl.DistributionParameters); 
Mu = Params(2*(1:2)-1,1:2); % Extract the means
Sigma = zeros(2,2,2);
for j = 1:2
    Sigma(:,:,j) = diag(Params(2*j,:)).^2; % Create diagonal covariance matrix
    xlim = Mu(j,1) + 4*[-1 1]*sqrt(Sigma(1,1,j));
    ylim = Mu(j,2) + 4*[-1 1]*sqrt(Sigma(2,2,j));
    f = @(x1,x2)reshape(mvnpdf([x1(:),x2(:)],Mu(j,:),Sigma(:,:,j)),size(x1));
    fcontour(f,[xlim ylim]) % Draw contours for the multivariate normal distributions 
end
% h.XLim = cxlim;
% h.YLim = cylim;
h.XLim = [0.39571,1.4542];
h.YLim = [0.1500, 0.7500];
title(' Expectation Maximization')
ylabel('Eccentricity')
xlabel('Hue')
legend('apple','orange')
colorbar();
hold off
