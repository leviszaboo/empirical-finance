%% Stock Data Pepsico Inc. (2000-2022)
% For Empirical Finance 3.2
% Assignment 1
% By team 25

addpath('./hist_stock_data');

f=0;

%% Question 1 (a)
% Load the data using hist_stock_data.m
price = hist_stock_data('01012000', '31122022', 'PEP');

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

%% Question 2 (a)

% (a) Fit an AR(1)-ARCH(2) model and display the output
MdlAR1ARCH2 = arima('ARLags', 1, 'Variance', garch(1, 2));

% Estimate the AR(1)-ARCH(2) model
EstMdlAR1ARCH2 = estimate(MdlAR1ARCH2, simple_y);

% Display the estimated coefficients
disp('AR(1)-ARCH(2) Coefficients:');
disp(EstMdlAR1ARCH2.AR);
disp(EstMdlAR1ARCH2.Variance.ARCH);

%% Question 2 (b) Compute standardized residuals and perform Jarque-Bera test
[resAR1ARCH2, varAR1ARCH2, logLAR1ARCH2] = infer(EstMdlAR1ARCH2, simple_y);
stdresAR1ARCH2 = resAR1ARCH2 ./ sqrt(varAR1ARCH2);

% Jarque-Bera test
[h_JB_AR1ARCH2, p_JB_AR1ARCH2, jbstat_JB_AR1ARCH2, critval_JB_AR1ARCH2] = jbtest(stdresAR1ARCH2);

% Display the Jarque-Bera test results
disp('Jarque-Bera Test for AR(1)-ARCH(2) Residuals:');
disp(['Hypothesis Test (H0: Residuals are normally distributed): ', num2str(h_JB_AR1ARCH2)]);
disp(['P-value: ', num2str(p_JB_AR1ARCH2)]);
disp(['J-B statistic: ', num2str(jbstat_JB_AR1ARCH2)]);
disp(['Critical values at the 5% significance level: ', num2str(critval_JB_AR1ARCH2)]);

%% Question 2 (c) Fit an AR(2)-GARCH(1,2) model and display the output
MdlAR2GARCH12 = arima('ARLags', 2, 'Variance', garch(1, 2));

% Estimate the AR(2)-GARCH(1,2) model
EstMdlAR2GARCH12 = estimate(MdlAR2GARCH12, simple_y);

% Display the estimated coefficients
disp('AR(2)-GARCH(1,2) Coefficients:');
disp(EstMdlAR2GARCH12.AR);
disp(EstMdlAR2GARCH12.Variance.ARCH);

%% Question 2 (d) Compute standardized residuals and perform Jarque-Bera test for AR(2)-GARCH(1,2)
[resAR2GARCH12, varAR2GARCH12, logLAR2GARCH12] = infer(EstMdlAR2GARCH12, simple_y);
stdresAR2GARCH12 = resAR2GARCH12 ./ sqrt(varAR2GARCH12);

% Jarque-Bera test
[h_JB_AR2GARCH12, p_JB_AR2GARCH12, jbstat_JB_AR2GARCH12, critval_JB_AR2GARCH12] = jbtest(stdresAR2GARCH12);

% Display the Jarque-Bera test results
disp('Jarque-Bera Test for AR(2)-GARCH(1,2) Residuals:');
disp(['Hypothesis Test (H0: Residuals are normally distributed): ', num2str(h_JB_AR2GARCH12)]);
disp(['P-value: ', num2str(p_JB_AR2GARCH12)]);
disp(['J-B statistic: ', num2str(jbstat_JB_AR2GARCH12)]);
disp(['Critical values at the 5% significance level: ', num2str(critval_JB_AR2GARCH12)]);

%% (e) Likelihood Ratio Test
% Loglikelihood value for AR(1)-ARCH(2) model
[~, ~, logL_AR1ARCH2] = infer(EstMdlAR1ARCH2, simple_y);

% Loglikelihood value for AR(2)-GARCH(1,2) model
[~, ~, logL_AR2GARCH12] = infer(EstMdlAR2GARCH12, simple_y);

% Perform Likelihood Ratio Test
dof = 1;  % Difference in the number of parameters between models
[h1,pValue1,stat1,cValue1] = lratiotest(logL_AR1ARCH2, logL_AR2GARCH12, dof)

%% Question 3 (a)

% Fit AR(2)-GARCH(1,2) model with t-distributed errors
MdlAR2G12 = arima('ARLags', [1, 2], 'Variance', garch('ARCHLags', 1, 'GARCHLags', [1, 2], 'Distribution', 't'));
EstMdlAR2G12 = estimate(MdlAR2G12, simple_y);

% Obtain standardized residuals
[resAR2G12, varAR2G12, logLAR2G12] = infer(EstMdlAR2G12, simple_y);
stdresAR2G12 = resAR2G12 ./ sqrt(varAR2G12);

% Plot return residuals
f=f+1; 
figure(f)
hold on
plot(price.Date(2:end), stdresAR2G12)
ylabel('Standardized Residuals')
legend('GARCH(1,2) Standardized Residuals')
xlabel('Date');
hold off 

%% Question 3 (b)

% Fit AR(2)-GARCH(1,2) model with standard normally distributed errors
MdlNormal = arima('ARLags', [1, 2], 'Variance', garch('ARCHLags', 1, 'GARCHLags', [1, 2], 'Distribution', 'Gaussian'));
EstMdlNormal = estimate(MdlNormal, simple_y);


% Obtain loglikelihood
[~, ~, logLNormal] = infer(EstMdlNormal, simple_y);

% Perform Likelihood Ratio Test
logLU = logLAR2G12; % Log-likelihood of unrestricted model
logLR = logLNormal; % Log-likelihood of restricted model
dof = 1;  % Difference in number of parameters
[h2,pValue2,stat2,cValue2] = lratiotest(logLU,logLR,dof)

%% Question 4 (a)

% Although in this section, the for loops, storing a single date in an array 
% and the array of sigmas and VAR values are not convenient for one forecast,
% it makes the code easily reusable if we include additional dates.

% Last trading day of 2022
str=["2022-12-30"];

% Conversion into datetime and find matching date indexes
date_arr = datetime(str,'InputFormat','yyyy-MM-dd', 'Format','dd/MM/yyyy');
dates=datetime(price.Date,'Format','dd/MM/yyyy');

[logical_index, date_arr_index]=ismember(date_arr, dates);

date_arr_index=date_arr_index-1;

% Event windows
W_E=[500 700 1000];

% Calculate daily VaR of Pepsico at 5% for 100€ investment
p=0.05;
val=100;
VAR=zeros(length(W_E),length(date_arr)); 

for i=1:length(date_arr_index)
    for j=1:length(W_E)
        date_arr_index(i)-W_E(j)+1;
        data=simple_y(date_arr_index(i)-W_E(j)+1:date_arr_index(i)); 
        data=sort(data); 
        VAR(j,i)= -data(W_E(j)*p)*val;
    end 
end 

VAR

%% Question 4 (b)

sigma = zeros(length(date_arr),1);

for i=1:length(date_arr) 

    % Estimating the model
    EstMdlG12 = estimate(MdlAR2G12, simple_y(1:date_arr_index(i)));
    [resG12,varG12,logLG12] = infer(EstMdlG12, simple_y(1:date_arr_index(i)));

    const = EstMdlG12.Variance.Constant
    alpha = cell2mat(EstMdlG12.Variance.ARCH)
    beta = cell2mat(EstMdlG12.Variance.GARCH)

    % Forecast parameters
    epsilon = resG12(date_arr_index(i))
    variance = varG12(date_arr_index(i))
    variance_lag1 = varG12(date_arr_index(i) - 1)

    % Computing variance forecast

    sigma(i) = sqrt(const + alpha*epsilon^2+beta(1)*variance+beta(2)*variance_lag1);
end 

VAR_Garch=zeros(length(date_arr),1);

for i=1:length(date_arr_index)
    VAR_Garch(i) = -sigma(i)*norminv(p)*val;
end 

VAR_Garch

%% Question 4 (c)

% (I) Simple ES

ES_Simple = zeros(length(W_E),length(date_arr));

% Calculating the average of all observations less than or equal to p ∗ W * 100 

for i=1:length(date_arr_index)
    for j=1:length(W_E)
        data=simple_y(date_arr_index(i)-W_E(j)+1:date_arr_index(i)); 
        data=sort(data); 
        ES_Simple(j,i)= -mean(data(1:W_E(j)*p))*val;
    end 
end 

ES_Simple

% (I) GARCH(1,2)-based ES

ES_Garch= zeros(length(date_arr),1); 

for i=1:length(date_arr_index)
    ES_Garch(i) = sigma(i)*normpdf(norminv(p))*val/p;
end 

ES_Garch





