%Chapter 6 - HAYASHI - EXERCISE 1: a,b,c,d,e,f,g
%Yen and FED FRED has been included in this folder
%John Daniel Paletto
%Last modified 12/10/15





%!@#$%^&*()PLEASE REFRAIN FROM CLICKING "RUN", PLEASE USE "RUN SECTION"!@#$%^&*()

%%%%% THIS WILL RUN PART BY PART and avoid DEATH-BY-POP-UPS





clear       %delete/clear memory
clc         %clear output screen
close all   %close e.g. figures

load('yen.mat'); %Jap Yen
%Notes from [currency].mat - Column 1,2,3,4: Date, Ask price S(t) in spot, 30-day
%forward F(t), bid price in delivery date for forward contract in spot
%ALL IN UNITS OF FOREIGN CURRENCY

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identify the week where forward premium is largest
%For that week find a 1-mo measure of the interst rate (US and Foreign)
%Verify Forward premium

%initialize
n=length(yen);
st=log(yen(:,2));
ft=log(yen(:,3));
s30t=log(yen(:,4));

%Calculate Forward Premium
fwd_premium=ft-st;
%annualize
fwd_premium=fwd_premium*12;

%From Federal Reserve Bank of St. Louis:
%INTGSTJPM193N Interest Rates, Government Securities, Government Bonds for Japan
%INTGSTUSM193N Interest Rates, Government Securities, Treasury Bills for United States
load('fred.mat');

%Calculate difference in Interest rates
dif=fredgraph.INTGSTJPM193N-fredgraph.INTGSTUSM193N;

%Find max
fwdmax=max(abs(fwd_premium));
[row,column] = find(fwd_premium==fwdmax);
if length(row)==0
    [row,column] = find(fwd_premium==-fwdmax);
    fwdmax=-fwdmax;
end;

difmax=max(abs(dif));
[difrow,difcolumn] = find(dif==difmax);

if length(difrow)==0
    [difrow,difcolumn] = find(dif==-difmax);
    difmax=-difmax;
end;

%OUTPUT

%Plots
startDate = datenum('01-01-1975');
endDate = datenum('12-01-1989');
xData = linspace(startDate,endDate,180);
xxData = linspace(startDate,endDate,778);

figure
plot(xData, dif);
hold on;
plot(xData(difrow),difmax,'*');
hold on;
plot(xxData,100*fwd_premium);
hold on;
plot(xxData(row),100*fwdmax,'*');
datetick('x','yyyy','keeplimits')

title('Difference in interest rates vs Forward Premium')
legend('i*-i','Difference in interest Max','Yen Foward Premium','FWD Premium Max')
ylabel('%')
xlabel('t')
datetick('x','yyyy','keeplimits')
hold off;

figure
plot(xData, fredgraph.INTGSTUSM193N);
hold on;
plot(xData, fredgraph.INTGSTJPM193N);
legend('US 1mo Gov Securities T-Bills','JP 1mo Gov Securities T-Bills')
datetick('x','yyyy','keeplimits')
title('Interest Rates')
ylabel('%')
xlabel('t')
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw sample correlogram of eps with 40 lags
%Significance of four lags?
%Check mean, and verify the assumption: Mean zero
%Subtract observed mean

%Calculate error
e = (s30t - ft);
%Annualize
e = 12*100*e;

%Autocovariance and auto correlation
rho_e=autocorrelnm(e,40);
gamma_e=autocov(e,40);

%standard error
std_e=1/sqrt(n);

%PLOT of forcast error
figure
plot(xxData,e)
datetick('x','yyyy','keeplimits')
title('Forecast Error')
ylabel('e')
xlabel('t')
hold off;

%PLOT correlogram
x=1:40;
zero=zeros(1,40);
up_errorband=zeros(1,40)+std_e*2;
low_errorband=zeros(1,40)-std_e*2;

figure
plot(x,rho_e,'b-o',x,zero,'g--')
hold on;
plot(x,up_errorband,'r--',x,low_errorband,'r--')
title('Correlogram of e - 40 lags - No mean subtracted')
legend('','','2 standard errors')
ylabel('rho')
xlabel('lags')
hold off;

%OUTPUT
MEAN_OF_e=mean(e)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Verify if currency is a random walk with drift
%Correlogram S(t+1)-S(t) with 40 lags
%Box-Ljung Q hypothesis for S(t) is a random walk with drift

%Calculate s(t+1) - s(t) for 40 observations
s1=st(1:end-1)-st(2:end);

%Autocovariance and auto correlation
rho_s1=autocorrelnm(s1,40);
gamma_s1=autocov(s1,40);

%Correlogram
x=1:40;
zero=zeros(1,40);
up_errorband=zeros(1,40)+std_e*2;
low_errorband=zeros(1,40)-std_e*2;

figure
plot(x,rho_s1,'b-o',x,zero,'g--')
hold on;
plot(x,up_errorband,'r--',x,low_errorband,'r--')
legend('','','2 standard errors')
title('Correlogram of s(t+1)-s(t) - 40 lags')
ylabel('rho')
xlabel('lags')
hold off;

%LJUNG Q test
q_s1 = LjungQ( s1, 40 );

for lag=1:40
critical_value_Ljung = icdf('chi2',0.95,lag);
if q_s1(lag)<critical_value_Ljung
    lag
    'based on Ljung´s test for conditional heteroskedasticity, we fail to reject the assumption of The data are independently distributed'
    
else
    lag
    'based on Ljungs´s test for conditional heteroskedasticity, we reject the assumption of independently distributed data; they exhibit serial correlation'
        
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%DICKEYFULLER & LAG LENGTH%%%%%%%%%%%%%%%%%%%%%%%
%Test S(t) for a unit root using Augmented DF (with and without trend/intrcpt)
%t-testing
%AIC
%BIC
%Verify choices, choose most appropriate

P=st;
T=length(st);

%plot the series to decide which type of test to use
figure,
plot(P)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TYPE I TEST: No intercept/no trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=P;                                            %select the series that needs to be tested
maxlag=floor(12*((T/100)^(1/4)));               %p_max based on 'rule of thumb'
t_P_1_xi=zeros(maxlag,1);                       %general-to-specific test statistic (check the last xi in each iteration for significance)
BIC_P_1=zeros(maxlag,1);                        %Bayesian-Schwartz Information Criterion
AIC_P_1=zeros(maxlag,1);                        %Akaike Information Criterion
rho_P_1=zeros(maxlag,1);                        %estimate for rho
t_P_1_rho=zeros(maxlag,1);                      %ADF t-statistic 
indicator_t=zeros(maxlag,1);                      %indicates whether the last xi estimate is significant
for i=maxlag:-1:1 %go from "big" model to "small" model
    indep=[];
    dep=y(i+2:end,:); %generate left-hand side variable
    for j=1:i
        indep=[indep y(i+2-j:end-j,:)-y(i+1-j:end-j-1,:)]; %generate set of regressors
    end;
    nn=length(indep); %notice: you lose observations as you increase the # of lags
    indep=[y(i+1:end-1,:) indep];
    
    coeff=indep\dep; %OLS regression
    rho_P_1(i,1)=coeff(1,1); %store the rho coefficient
    resid=dep-indep*coeff; %OLS residuals
    
    %estimate for the variance covariance matrix of errors
    s2e=(1/(nn-1-i))*(resid'*resid); 
    
    %Avar
    help=((indep'*indep)^(-1));
    SE=diag((s2e*help).^(1/2));
    clear help s2e
    
    %test for insigificance of the last autoregressive coefficient 
    t_P_1_xi(i,1)=(coeff(end,1)-0)/SE(end,1);
    %ADF t statistic
    t_P_1_rho(i,1)=(coeff(1,1)-1)/SE(1,1);
    clear SE
    
    %Critical value for the t test for general-to-specific approach
    critical_value_t=icdf('Normal',0.975,0,1);
    
    if abs(t_P_1_xi(i,1))>=critical_value_t
        indicator_t(i,1)=1;
    end;
    
    %information criteria (NOTE THAT WE ARE EFFECTIVELY ESTIMATING AN
    %AR(p+1) on the level of y!!!
    %we use the "T-maxlag" fixed sample size
    BIC_P_1(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*((log((T-maxlag-1)))/(T-maxlag-1));
    AIC_P_1(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*(2/(T-maxlag-1));
end;
clear i indep dep j nn coeff resid Sigma_e help Avar coeff2 mm R r critical_value_t

%for p=0
indep=[];
dep=y(2:end,1);
indep=y(1:end-1,1);

coeff=indep\dep; %OLS regression
resid=dep-indep*coeff;
s2e=(1/(T-1-1))*(resid'*resid); 
help=((indep'*indep)^(-1));
SE=diag((s2e*help).^(1/2));
clear help s2e

%all these are sorted such that the first element corresponds to the
%largest model and the last element to the smallest model
rho_P_1=[coeff(1,1); rho_P_1];
t_P_1_rho=[(coeff(1,1)-1)/SE(1,1); t_P_1_rho];
BIC_P_1=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*((log((T-maxlag-1)))/(T-maxlag-1)); BIC_P_1];
AIC_P_1=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*(2/(T-maxlag-1)); AIC_P_1];
clear indep dep coeff resid y SE

[min_AIC,indicator_AIC]=min(AIC_P_1);
[min_BIC,indicator_BIC]=min(BIC_P_1);

fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'\n\n                      Result from a type I (A)DF Test                          \n');
    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'The DF-t test statistic is equal to:  \t %f       \n', t_P_1_rho(1,1));
    fprintf(1,'The OLS estimate for rho is equal to: \t %f \n', rho_P_1(1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The general-to-specific lag selection approach selects  \t %8.0f       \t lags \n', max(find(indicator_t==1)));
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f       \n', t_P_1_rho(max(find(indicator_t==1))+1,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_1(max(find(indicator_t==1))+1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The BIC information criterion selects  %8.0f       \t lags \n', indicator_BIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to: \t %f      \n', t_P_1_rho(indicator_BIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_1(indicator_BIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The AIC information criterion selects  %8.0f       \t lags \n', indicator_AIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f      \n', t_P_1_rho(indicator_AIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_1(indicator_AIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
clear indicator_t min_AIC indicator_AIC min_BIC indicator_BIC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TYPE II TEST: Intercept/no trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=P;                                            %select the series that needs to be tested
maxlag=floor(12*((T/100)^(1/4)));               %p_max based on 'rule of thumb'
t_P_2_xi=zeros(maxlag,1);                       %general-to-specific test statistic (check the last xi in each iteration for significance)
BIC_P_2=zeros(maxlag,1);                        %Bayesian-Schwartz Information Criterion
AIC_P_2=zeros(maxlag,1);                        %Akaike Information Criterion
rho_P_2=zeros(maxlag,1);                        %estimate for rho
t_P_2_rho=zeros(maxlag,1);                      %ADF t-statistic 
indicator_t=zeros(maxlag,1);                      %indicates whether the last xi estimate is significant
for i=maxlag:-1:1 %go from "big" model to "small" model
    indep=[];
    dep=y(i+2:end,:); %generate left-hand side variable
    for j=1:i
        indep=[indep y(i+2-j:end-j,:)-y(i+1-j:end-j-1,:)]; %generate set of regressors
    end;
    nn=length(indep); %notice: you lose observations as you increase the # of lags
    indep=[ones(nn,1) y(i+1:end-1,:) indep];
    
    coeff=indep\dep; %OLS regression
    rho_P_2(i,1)=coeff(2,1); %store the rho coefficient
    resid=dep-indep*coeff; %OLS residuals
    
    %estimate for the variance covariance matrix of errors
    s2e=(1/(nn-2-i))*(resid'*resid); 
    
    %Avar
    help=((indep'*indep)^(-1));
    SE=diag((s2e*help).^(1/2));
    clear help s2e
    
    %test for insigificance of the last autoregressive coefficient 
    t_P_2_xi(i,1)=(coeff(end,1)-0)/SE(end,1);
    %ADF t statistic
    t_P_2_rho(i,1)=(coeff(2,1)-1)/SE(2,1);
    clear SE
    
    %Critical value for the t test for general-to-specific approach
    critical_value_t=icdf('Normal',0.975,0,1);
    
    if abs(t_P_2_xi(i,1))>=critical_value_t
        indicator_t(i,1)=1;
    end;
    
    %information criteria (NOTE THAT WE ARE EFFECTIVELY ESTIMATING AN
    %AR(p+1) on the level of y!!!
    %we use the "T-maxlag" fixed sample size
    BIC_P_2(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*((log((T-maxlag-1)))/(T-maxlag-1));
    AIC_P_2(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*(2/(T-maxlag-1));
end;
clear i indep dep j nn coeff resid Sigma_e help Avar coeff2 mm R r critical_value_t

%for p=0
indep=[];
dep=y(2:end,1);
indep=[ones(T-1,1) y(1:end-1,1)];

coeff=indep\dep; %OLS regression
resid=dep-indep*coeff;
s2e=(1/(T-2-1))*(resid'*resid); 
help=((indep'*indep)^(-1));
SE=diag((s2e*help).^(1/2));
clear help s2e

%all these are sorted such that the first element corresponds to the
%largest model and the last element to the smallest model
rho_P_2=[coeff(2,1); rho_P_2];
t_P_2_rho=[(coeff(2,1)-1)/SE(2,1); t_P_2_rho];
BIC_P_2=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*((log((T-maxlag-1)))/(T-maxlag-1)); BIC_P_2];
AIC_P_2=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*(2/(T-maxlag-1)); AIC_P_2];
clear indep dep coeff resid y SE

[min_AIC,indicator_AIC]=min(AIC_P_2);
[min_BIC,indicator_BIC]=min(BIC_P_2);

fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'\n\n                      Result from a type II (A)DF Test                          \n');
    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'The DF-t test statistic is equal to:  \t %f       \n', t_P_2_rho(1,1));
    fprintf(1,'The OLS estimate for rho is equal to: \t %f \n', rho_P_2(1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The general-to-specific lag selection approach selects  \t %8.0f       \t lags \n', max(find(indicator_t==1)));
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f       \n', t_P_2_rho(max(find(indicator_t==1))+1,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_2(max(find(indicator_t==1))+1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The BIC information criterion selects  %8.0f       \t lags \n', indicator_BIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to: \t %f      \n', t_P_2_rho(indicator_BIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_2(indicator_BIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The AIC information criterion selects  %8.0f       \t lags \n', indicator_AIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f      \n', t_P_2_rho(indicator_AIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_2(indicator_AIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
clear indicator_t min_AIC indicator_AIC min_BIC indicator_BIC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TYPE III TEST: Intercept/trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=P;                                            %select the series that needs to be tested
maxlag=floor(12*((T/100)^(1/4)));               %p_max based on 'rule of thumb'
t_P_3_xi=zeros(maxlag,1);                       %general-to-specific test statistic (check the last xi in each iteration for significance)
BIC_P_3=zeros(maxlag,1);                        %Bayesian-Schwartz Information Criterion
AIC_P_3=zeros(maxlag,1);                        %Akaike Information Criterion
rho_P_3=zeros(maxlag,1);                        %estimate for rho
t_P_3_rho=zeros(maxlag,1);                      %ADF t-statistic 
indicator_t=zeros(maxlag,1);                      %indicates whether the last xi estimate is significant
for i=maxlag:-1:1 %go from "big" model to "small" model
    indep=[];
    dep=y(i+2:end,:); %generate left-hand side variable
    for j=1:i
        indep=[indep y(i+2-j:end-j,:)-y(i+1-j:end-j-1,:)]; %generate set of regressors
    end;
    nn=length(indep); %notice: you lose observations as you increase the # of lags
    indep=[ones(nn,1) [1:1:nn]' y(i+1:end-1,:) indep];  %generate set of regressors
    
    coeff=indep\dep; %OLS regression
    rho_P_3(i,1)=coeff(3,1); %store the rho coefficient
    resid=dep-indep*coeff; %OLS residuals
    
    %estimate for the variance covariance matrix of errors
    s2e=(1/(nn-3-i))*(resid'*resid); 
    
    %Avar
    help=((indep'*indep)^(-1));
    SE=diag((s2e*help).^(1/2));
    clear help s2e
    
    %test for insigificance of the last autoregressive coefficient 
    t_P_3_xi(i,1)=(coeff(end,1)-0)/SE(end,1);
    %ADF t statistic
    t_P_3_rho(i,1)=(coeff(3,1)-1)/SE(3,1);
    clear SE
    
    %Critical value for the t test for general-to-specific approach
    critical_value_t=icdf('Normal',0.975,0,1);
    
    if abs(t_P_3_xi(i,1))>=critical_value_t
        indicator_t(i,1)=1;
    end;
    
    %information criteria (NOTE THAT WE ARE EFFECTIVELY ESTIMATING AN
    %AR(p+1) on the level of y!!!
    %we use the "T-maxlag" fixed sample size
    BIC_P_3(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*((log((T-maxlag-1)))/(T-maxlag-1));
    AIC_P_3(i,1)=log((1/(T-maxlag-1))*(resid(maxlag-i+1:end,1)'*resid(maxlag-i+1:end,1)))+(i+1+1)*(2/(T-maxlag-1));
end;
clear i indep dep j nn coeff resid Sigma_e help Avar coeff2 mm R r critical_value_t

%for p=0
indep=[];
dep=y(2:end,1);
indep=[ones(T-1,1) [1:1:T-1]' y(1:end-1,1)];

coeff=indep\dep; %OLS regression
resid=dep-indep*coeff;
s2e=(1/(T-3-1))*(resid'*resid); 
help=((indep'*indep)^(-1));
SE=diag((s2e*help).^(1/2));
clear help s2e

%all these are sorted such that the first element corresponds to the
%largest model and the last element to the smallest model
rho_P_3=[coeff(3,1); rho_P_3];
t_P_3_rho=[(coeff(3,1)-1)/SE(3,1); t_P_3_rho];
BIC_P_3=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*((log((T-maxlag-1)))/(T-maxlag-1)); BIC_P_3];
AIC_P_3=[log((1/(T-maxlag-1))*(resid(maxlag+1:end,1)'*resid(maxlag+1:end,1)))+(0+1+1)*(2/(T-maxlag-1)); AIC_P_3];
clear indep dep coeff resid y SE

[min_AIC,indicator_AIC]=min(AIC_P_3);
[min_BIC,indicator_BIC]=min(BIC_P_3);

fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'\n\n                      Result from a type III (A)DF Test                          \n');
    fprintf(1,'\n\n-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'-----------------------------------------------------------------------------------------------------\n\n');
    fprintf(1,'The DF-t test statistic is equal to:  \t %f       \n', t_P_3_rho(1,1));
    fprintf(1,'The OLS estimate for rho is equal to: \t %f \n', rho_P_3(1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The general-to-specific lag selection approach selects  \t %8.0f       \t lags \n', max(find(indicator_t==1)));
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f       \n', t_P_3_rho(max(find(indicator_t==1))+1,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_3(max(find(indicator_t==1))+1,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The BIC information criterion selects  %8.0f       \t lags \n', indicator_BIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to: \t %f      \n', t_P_3_rho(indicator_BIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_3(indicator_BIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n');
    fprintf(1,'The AIC information criterion selects  %8.0f       \t lags \n', indicator_AIC-1);
    fprintf(1,'The corresponding ADF-t test statistic is equal to:  \t %f      \n', t_P_3_rho(indicator_AIC,1));
    fprintf(1,'The corresponding OLS estimate for rho is equal to: \t %f \n', rho_P_3(indicator_AIC,1));
    fprintf(1,  '-----------------------------------------------------------------------------------------------------\n\n');
clear indicator_t min_AIC indicator_AIC min_BIC indicator_BIC


%%%%%%%%%%%%%%%%%%%%%%%%%%%END DICKEY FULLER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unconditional test
%Replicate Table 6.1

%Exchange rates
actual_rate=(s30t-st)*12*100;
expected_rate=(ft-st)*12*100;
unexpected=e;

%Means
mean_actual=mean(actual_rate);
mean_expected=mean(expected_rate);
mean_e=MEAN_OF_e;


%Standard deviations
stdd_actual=std(actual_rate);
stdd_expected=std(expected_rate);
stdd_e=std(e);

%Standard error from proposition 6.10
std_e61=sqrt((gamma_e(1)+sum(gamma_e(2:5))*2))/n;

%Table 6.1
actual=[mean_actual,stdd_actual,0];
expected=[mean_expected,stdd_expected,0];
unex=[mean_e, stdd_e, std_e61];

Table61(:,1)=actual;
Table61(:,2)=expected;
Table61(:,3)=unex;
cnames = {'s30 - s','f - s','Difference'};
rnames = {'Mean','Std Deviation','Std Error'};
set(figure,'name', 'Table 6.1, Y/$ - Means and Standard Deviations','numbertitle','off');
uitable('Data',Table61,'ColumnName',cnames,'RowName',rnames,'Position',[20 20 335 335])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Regression test with truncated kernel
%Replicate Table 6.2

X(:,1)=ones(1,778);
X(:,2)=expected_rate;
[T,K]=size(X);

y=actual_rate;

%SHORT SOLUTION: 
mdl = fitlm(expected_rate,y,'linear')
[EstCovtr,LSSe,coeff]=hac(mdl,'type','HAC','weights','TR','bandwidth',4,'smallT',false)
[h,p,s,cv]=waldtest([coeff(1);coeff(2)-1],[1 0; 0 1],EstCovtr,.01)
%%

%LONG SOLUTION
delta_hat_GMM=X\y;
epsilon_hat=y-X*delta_hat_GMM;
clear delta_hat_OLS

%define g (i.e. the multiplication of residuals*regressors - or
%residuals*instruments)
g_hat=X.*repmat(epsilon_hat,[1,K]);
clear epsilon_hat

hat_Gamma_j=NaN(K,K,2*(T-1)+1);
for j=0:1:(T-1) 
    help=0; %auxiliary variable
    for t=j+1:1:T %summation index
        help=help+g_hat(t,:)'*g_hat(t-j,:);
    end;
    clear t
    hat_Gamma_j(:,:,T+j)=(1/T)*help; %compute Gamma_j for lags 0 to T-1
    if j>0
        hat_Gamma_j(:,:,T-j)=(reshape(hat_Gamma_j(:,:,T+j),[K,K]))'; %compute remaining Gamma_j's for lags -1 to -(T-1)
    end;
    clear help
end;
clear j 

%Given
q=4;

%For the sake of exhibition, take the same bandwidth/window size as above

Omega_hat_Truncated=zeros(K,K);
kernel_Truncated=NaN(2*(T-1)+1,1);
for j=-(T-1):1:(T-1)
    x=(j/q); %this is the kernel arguement j/q(T)
    if abs(x)<=1
        kernel_Truncated(T+j,1)=1;
    else
        kernel_Truncated(T+j,1)=0;
    end;
    Omega_hat_Truncated=Omega_hat_Truncated+kernel_Truncated(T+j,1)*(reshape(hat_Gamma_j(:,:,T+j),[K,K]));
    clear x
end;
clear j 

Sxx = X'*X/(n);
sxy = X'*y/(n);
delta_OLS=(Sxx^(-1))*sxy;
e_hat=y-X*delta_OLS;

Avar_GMM_robust=(Sxx'*(Omega_hat_Truncated^(-1))*Sxx)^(-1);
SE_GMM_robust=diag(((1/T)*Avar_GMM_robust).^(1/2));

R=[1, 0];
r=1;


R2=1-(e_hat'*e_hat)/((y-mean(y))'*(y-mean(y)));
SE_R = ((e_hat-mean(e_hat))'*(e_hat-mean(e_hat))/(n-2))^(1/2);

figure,
plot(expected_rate,actual_rate,'o')
hold on;
plot(expected_rate,delta_OLS(1)+delta_OLS(2)*expected_rate)
title('Fig. 6.5: Regression of Actual against Expected Rates')
legend('Scatter','Bo + B1(f-s)')
ylabel('s30-s')
xlabel('f-s')
hold off;

R_tbl(1,1:2)=delta_OLS;
R_tbl(2,1:2)=[SE_GMM_robust(1),SE_GMM_robust(2)];
R_tbl(3,3:5)=[R2,mean(y),s];

freg=figure('Position', [150 150 600 255]);
set(freg, 'name', 'TABLE 6.2: Regression Tests of Market Efficiency: 1975-1989',...
    'numbertitle','off');

r2names={'Coefficients','Std Error','Statistics'};
c2names={,'B0','B1','R^2','Mean of y', 'Wald-stat'};
Regression_table = uitable('Data',R_tbl,...
    'RowName',r2names,'ColumnName',c2names,'Tag',...
    'Regression Tests of Market Efficiency: 1975-1989',...
    'Parent', freg,'Position',[40 40 550 155]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End Part E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART F%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bartlett kernel based estimator of S for regression
%Newey-West data-dependent automatic bandwidth selection
%Assume lag length is 12 for YEN

%SHORT SOLUTION:
[EstCovbt,LSSe,coeff]=hac(mdl,'type','HAC','weights','BT','bandwidth',13,'smallT',false)
[h,p,s,cv]=waldtest([coeff(1);coeff(2)-1],[1 0; 0 1],EstCovbt,.01)
%%

%LONG SOLUTION
%GIVEN
q=12+1;

Omega_hat_Bartlett=zeros(K,K);
kernel_Bartlett=NaN(2*(T-1)+1,1);
for j=-(T-1):1:(T-1)
    x=(j/q); %this is the kernel arguement j/q(T)
    if abs(x)<=1
        kernel_Bartlett(T+j,1)=1-abs(x);
    else
        kernel_Bartlett(T+j,1)=0;
    end;
    Omega_hat_Bartlett=Omega_hat_Bartlett+kernel_Bartlett(T+j,1)*(reshape(hat_Gamma_j(:,:,T+j),[K,K]));
    clear x
end;
clear j 


Avar_GMM_bartlett_robust=(Sxx'*(Omega_hat_Bartlett^(-1))*Sxx)^(-1);
SE_GMM_bartlett_robust=diag(((1/T)*Avar_GMM_bartlett_robust).^(1/2));

R_tbl(1,1:2)=delta_OLS;
R_tbl(2,1:2)=[SE_GMM_bartlett_robust(1),SE_GMM_bartlett_robust(2)];
R_tbl(3,3:5)=[R2,mean(y),s];

freg=figure('Position', [150 150 600 255]);
set(freg, 'name', '(Bartlett Kernel)Regression Tests of Market Efficiency: 1975-1989',...
    'numbertitle','off');

r2names={'Coefficients','Std Error','Statistics'};
c2names={,'B0','B1','R^2','Mean of y', 'Wald-stat'};
Regression_table = uitable('Data',R_tbl,...
    'RowName',r2names,'ColumnName',c2names,'Tag',...
    'BARTLETT Regression Tests of Market Efficiency: 1975-1989',...
    'Parent', freg,'Position',[40 40 550 155]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART F%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START PART G%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VARHAC ESTIMATOR FOR S
%Justify p

p=int16(n^(1/3));%maximal number of lags in all equations


Phi_hat_temp=zeros(p*K,K); %temporary matrix of coefficients
residual=NaN(T,K); %residuals of VAR estimation
for k=1:K
    indep=[];
    dep=g_hat(p+1:end,k); %generate left-hand side variable
    for j=1:p
        indep=[indep g_hat(p+1-j:end-j,:)]; %generate set of regressors
    end;
    nn=length(indep); %notice: you lose observations as you increase the # of lags
    indep=[ones(nn,1) indep];
    
    reg=indep\dep; %OLS regression
    Phi_hat_temp(1:p*K,k)=reg(2:end,:); %we don't store the estimate for the intercept
    residual(T-nn+1:T,k)=dep-indep*reg; %OLS residuals
	clear indep dep j nn reg 
end;
clear k

%now store the coefficients in the proper format (see slides)
Phi_hat=NaN(K,K,p);
for ii=1:p
    Phi_hat(:,:,ii)=Phi_hat_temp(ii*K-(K-1):ii*K,:)';
end;
clear ii Phi_hat_temp

%Now compute Sigma_epsilon_hat
Sigma_resid_hat=zeros(K,K);
for t=p+1:1:T
    Sigma_resid_hat=Sigma_resid_hat+residual(t,:)'*residual(t,:);
end;
Sigma_resid_hat=(1/T)*Sigma_resid_hat;
clear t residual

%Construct Omega_hat_VARHAC 
temp=zeros(K,K);
for ii=1:p
    temp=temp+Phi_hat(:,:,ii);
end;
clear ii
Omega_hat_VARHAC=((eye(K,K)-temp)^(-1))*Sigma_resid_hat*((eye(K,K)-temp)^(-1))';
clear temp

Avar_GMM_varhac_robust=(Sxx'*(Omega_hat_VARHAC^(-1))*Sxx)^(-1)
SE_GMM_varhac_robust=diag(((1/T)*Avar_GMM_varhac_robust).^(1/2))
EstCovv=Avar_GMM_varhac_robust/T;

[h,p,s,cv]=waldtest([coeff(1);coeff(2)-1],[1 0; 0 1],EstCovv,.01)

R_tbl(1,1:2)=delta_OLS;
R_tbl(2,1:2)=[SE_GMM_varhac_robust(1),SE_GMM_varhac_robust(2)];
R_tbl(3,3:5)=[R2,mean(y),s];

freg=figure('Position', [150 150 600 255]);
set(freg, 'name', '(VARHAC)Regression Tests of Market Efficiency: 1975-1989',...
    'numbertitle','off');

r2names={'Coefficients','Std Error','Statistics'};
c2names={,'B0','B1','R^2','Mean of y', 'Wald-stat'};
Regression_table = uitable('Data',R_tbl,...
    'RowName',r2names,'ColumnName',c2names,'Tag',...
    'VARHAC Regression Tests of Market Efficiency: 1975-1989',...
    'Parent', freg,'Position',[40 40 550 155]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END PART G%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%