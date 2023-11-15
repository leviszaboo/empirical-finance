%% Stock Data Pepsico Inc. (2000-2022)
% For Empirical Finance 3.2
% Assignment 1
% By team 25

addpath('./hist_stock_data');

f=0;

%% Question 1 (a)
% Load the data using hist_stock_data.m
price = hist_stock_data('01012002', '31122022', 'PEP');

% % Plot returns 
f=f+1; 
figure(f)
hold on 
plot(price.Date, price.Close)
xlabel('Year')
ylabel('Price')
title('PepsiCo, Inc.')
hold off

%% Question 1 (b)
% Simple returns
simple_y=(price.Close(2:end)-price.Close(1:end-1))./price.Close(1:end-1);
% Continuously compounded returns (using natural log)
comp_y=diff(log(price.Close));

% Sanity check: correlation using corr(x,y) function 
corr(simple_y,comp_y)

% Creating two separate subplots for the two types of returns
f=f+1; 
figure(f)
hold on 
subplot(2,1,1),  plot(datetime(price.Date(2:end)),simple_y)
legend('Simple Returns')
xlabel('Year')
ylabel('Return')
title('PepsiCo Returns over Time')
subplot(2,1,2), plot(datetime(price.Date(2:end)),comp_y, 'Color', [1, 0.5, 0])
legend('Continuously Compounded Returns')
xlabel('Year')
ylabel('Return')
hold off

%% Question 1 (c) 
% Autocorrolation graph for PepsiCO squared returns   
[acf2,lags2,bounds2]=autocorr(simple_y.^2,250);

upper2=ones(length(acf2)-1,1)*bounds2(1);
lower2=ones(length(acf2)-1,1)*bounds2(2);

f=f+1; 
figure(f)
hold on
plot(lags2(2:end),acf2(2:end))
plot(lags2(2:end),upper2, 'color','red' )
plot(lags2(2:end),lower2, 'color','red')
legend('ACF Coefficient', '95% Confidence Bound')
xlabel('Nr. of Lags')
hold off

%% Question 1 (d) 
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',20)
% Robustness checks
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',25)
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',30)
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',35)
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',60)
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',120)
[h,p,stat,cValue] = lbqtest(simple_y,'Lags',200)




%% Question 1 (e)
[muhat,sigmahat]=normfit(simple_y);
standard_simple_y = (simple_y - muhat)./sigmahat;
% The QQ-plot of the returns

f=f+1; 
figure(f)
hold on
qqplot(standard_simple_y) 
hold off

% Sanity check: Jarque-Bèra test
[h,p,jbstat,critval] = jbtest(simple_y)

%% Question 3 (a)

% Fit AR(2)-GARCH(1,2) model with t-distributed errors
MdlAR2G12 = arima('ARLags', [1, 2], 'Variance', garch('ARCHLags', 1, 'GARCHLags', 2, 'Distribution', 't'));
EstMdlAR2G12 = estimate(MdlAR2G12, simple_y);

% Obtain standardized residuals
[resAR2G12, varAR2G12, logLAR2G12] = infer(EstMdlAR2G12, simple_y);
stdresAR2G12 = resAR2G12 ./ sqrt(varAR2G12);

% Plot return residuals
f=f+1; 
figure(f)
hold on
plot(price.Date(2:end), stdresAR2G12)
legend('GARCH(1,2) Residuals')
xlabel('Date');
ylabel('Standardized Returns');
hold off 

%% Question 3 (b)

% Fit AR(2)-GARCH(1,2) model with standard normally distributed errors
MdlNormal = arima('ARLags', [1, 2], 'Variance', garch('ARCHLags', 1, 'GARCHLags', 2, 'Distribution', 'Gaussian'));
EstMdlNormal = estimate(MdlNormal, simple_y);


% Obtain loglikelihood
[~, ~, logLNormal] = infer(EstMdlNormal, simple_y);

% Perform Likelihood Ratio Test
logLU = logLAR2G12;
logLR = logLNormal;
dof = 1;  
[h,pValue,stat,cValue] = lratiotest(logLU,logLR,dof)