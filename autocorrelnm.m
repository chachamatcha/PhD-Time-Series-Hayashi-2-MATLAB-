function [ acrl ] = autocorrelnm(x , maxlag)
%Takes in 1xn vector and calculates lagged autocorellations with zero mean
n = length(x);

acv= zeros(maxlag+1,1);

for h=0: maxlag
   acv(h+1)= x(1:n-h)' * x(1+h:n);
end

for i=1:maxlag
    acrl(i,:)=acv(i+1)/acv(1);
end