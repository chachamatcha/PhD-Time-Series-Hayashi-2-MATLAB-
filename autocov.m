function [ acv ] = autocov(x , maxlag)
%Takes in 1xn vector and calculates lagged autocov
n = length(x);
x= x - mean(x);
acv= zeros(maxlag+1,1);

for h=0: maxlag
   acv(h+1)= x(1:n-h)' * x(1+h:n);
end