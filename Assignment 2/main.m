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
par = 2000;
r=0.026;
yield =[5.34 4.12 3.72 3.7 3.86 3.92];
T= length(yield)

% Price the bond
coupon= r*par;
cash_flow=coupon+zeros(1,T)
cash_flow(6)=cash_flow(6)+par
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
