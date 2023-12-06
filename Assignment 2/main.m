%% Stock Data Intel Corporation (2012-2022)
% For Empirical Finance 3.2
% Assignment 2
% By team 40 (Levi Szabo, Barry Nicolas en Renzo Vermeulen)

%% House Keeping 

% Housekeeping
clear all;
close all;
clc; 

% Figure counter
f=0;

%% Question 1 (a)
par = 1000;
r= 0.073;
yield =[5.34 4.12 3.72 3.7 3.86 3.92 3.95 3.97];
T= length(yield)

% Price the bond
coupon= r*par;
cash_flow=coupon+zeros(1,T)
cash_flow(8)=cash_flow(8)+par
P= sum(cash_flow./((1+yield./100).^(1:T)))

%% Question 1 (b)
% Set seed for randomization
randn('state', 1512) % set the seed
% Draw of random numbers to add to the yield
S = 10; % number of simulations
sigma = 1.25 % daily yield volatility
eps = randn(1,S)* sigma; % generate random yield changes in a loop
% Add shift to the yield
ysim=zeros(T,S);
for s = 1:10
     ysim(:,s) = yield + eps(s);
end
% Plot the different yields
f=f+1;
figure(f)
hold on
plot(yield, '-','color', 'b', 'LineWidth',2)
plot(ysim,'--','color',[0.5 0.5 0.5])
legend('True','Simulated')
xlabel('Time')
ylabel('Yield (%)')
daspect([1 1.5 1]);
ylim([3 10])
axis tight; 
hold off

%% Question 2 (a)
% Load the data using hist_stock_data.m
price = hist_stock_data('01012012',['31122022'],'INTC');
% Calculating simple returns
simple_y=(price.Close(2:end)-price.Close(1:end-1))./price.Close(1:end-1)
% Obtain sample statistics for simple returns
mean(simple_y)
std(simple_y)
min(simple_y)
max(simple_y)
% Obtain first and third quartiles for simple returns
% Step 1 - rank the data
y = sort(simple_y);
% Compute first quartile
Q(1) = median(y(find(y<median(y))));
display (Q(1))
% Compute third quartile
Q(3) = median(y(find(y>median(y))));
display (Q(3))

%% Question 2 (b) p=0.01

% Prepare the estimations by setting up the Matlab environment
% VaR(0.01) !
T = length(simple_y);          % number of observations for return y
WE = 1000;                     % estimation window length
p = 0.01;                      % VaR probability
value = 1;                     % portfolio value

VaR = NaN(T-WE, 1);             % matrix to hold VaR forecasts for 2 models

% Compute the VaR forecasts for using historic simulations 
index = p*WE;

for t=1:(T-WE)
 data= simple_y(t:WE-1+t); 
 data=sort(data);
 VaR(t,1) = -1*data(index)*value;
end 

% Plot returns and VaR
f=f+1; 
figure(f)
hold on 
plot(datetime(price.Date(WE+2:end)), simple_y(WE+1:end),'color',[0.5 0.5 0.5]); 
plot(datetime(price.Date(WE+2:end)), (-1)*VaR(:,1));
xlabel('Time');
ylabel('Returns')
legend('','VaR(0.01) using HS');
hold off

%% Question 2 (b) p=0.05
%% Prepare the estimations by setting up the Matlab environment
% VaR(0.05) !
T = length(simple_y);          % number of observations for return y
WE = 1000;                     % estimation window length
p = 0.05;                      % VaR probability
value = 1;                     % portfolio value

VaR = NaN(T-WE,1);             % matrix to hold VaR forecasts for 1 model

% Compute the VaR forecasts using historic simulations 
index = p*WE;

for t=1:(T-WE)
 data= simple_y(t:WE-1+t); 
 data=sort(data);
 VaR(t,1) = -1*data(index)*value;
end 
% Save the model outcomes in a mat file
save tmp.mat

clear all;
close all;
clc; 

load tmp.mat
% Plot returns and VaR
f=f+1; 
figure(f)
hold on 
plot(datetime(price.Date(WE+2:end)), simple_y(WE+1:end),'color',[0.5 0.5 0.5]); 
plot(datetime(price.Date(WE+2:end)), (-1)*VaR(:,1));
xlabel('Time');
ylabel('Returns')
legend('','VaR(0.05) using HS');
hold off
