function Forecasting(df)
clc;

[Spec] = myfrc.SelectBestModel();
ROrder = Spec.P;
MOrder = Spec.Q;
I = 1;
ARcoef = Spec.ARcoef;
MAcoef = Spec.MAcoef;
ConstTerm = Spec.ConstTerm;
disp(df.finalSeries);
strModel = sprintf('Best Model : \tP=%d, \tI=%d , \tQ=%d', ROrder, I, MOrder);
disp(strModel);

disp('AR coefficient');
disp(ARcoef);
disp('MA coefficient');
disp(MAcoef);
disp('Constant Term');
disp(ConstTerm);

if ROrder==0
    ARcoef = [];
end
if MOrder==0
    MAcoef = [];
end

numPeriods = 16;
model = arima(ROrder, df.nDiff, MOrder);
model.Constant = NaN; %ConstTerm;
model.AR{1:ROrder} = ARcoef;
model.MA{1:MOrder} = MAcoef;
model.Variance = var(df.finalSeries);
v2 = var(df.finalSeries');

%model = arima('Constant', ConstTerm, 'AR', ARcoef, 'MA', MAcoef, 'Variance', V);
%model = arima('Constant', ConstTerm, 'AR', ARcoef, 'MA', MAcoef, 'Variance', V, 'P', ROrder, 'D', df.nDiff, 'Q', MOrder);
%fit = estimate(model, Y);

fit = estimate(model, df.finalSeries);
[X, YMSE] = forecast(fit, numPeriods, 'Y0', df.finalSeries);

Newdata = [df.finalSeries; MeanForecast];
figure(10);
plot(Newdata);

disp('MeanForecast');
disp(MeanForecast);

ForecastedSeries = MeanForecast;
disp('MeanRMSE');
disp(MeanRMSE);

disp('SigmaForecast');
disp(SigmaForecast);

disp('SigmaTotal');
disp(SigmaTotal);

[NewSeries]= df.inverseNormalDiff(ForecastedSeries);

[s1,c1]=size(inputData);
figure(7)
plot(NewSeries)
xlim([0,230])
set(gca,'XTick', 0:10:210)
set(gca,'XTickLabel', 0:10:210)
title('Forecasted Passenger Demand Merter Station 1July 7am-12 pm')

for i=1:numPeriods
    ForecastedTail(i, c1)=NewSeries(s1+i, c1);
end

disp('ForecastedTail');
disp(ForecastedTail);

[l , k] = size(NewSeries);
%disp(NewSeries1);

set(handles.ForcVal, 'string', num2str(NewSeries(l,k)));
end


function [ARcoef, MAcoef, ConstTerm, ARcoefSEE, MAcoefSEE, ConstTermSEE, Residuals, IsStatInvert, IsConverged] = Estimate(Data, p, q)
%set initial parameters for AR 
ARcoef = [];
MAcoef = [];
ARcoefSEE = [];
MAcoefSEE = [];
for i=1:p
    ARcoef(i)=0.1;
end

%set initial parameters for MA
MAcoef = [];
for i=1:q
    MAcoef(i)=0.1;
end

%compute estimate for constant term
SumOfARCoef = sum(ARcoef);
ConstTerm = mean(Data)*(1-SumOfARCoef); 

%Model specification generation
ModelSpec = garchset('C', ConstTerm, 'R', p, 'M',q, 'AR', ARcoef, 'MA', MAcoef, 'Display','on');

%Coefficient Estimation
[Coefficients, Errors, LLF,Residuals, Sigma, Summery]= garchfit(ModelSpec, Data);
disp(Coefficients);

%extracting values
if p > 0
    ARcoef = Coefficients.AR; % garchget(Coefficients, 'AR');
end

if q > 0
    MAcoef = garchget(Coefficients, 'MA');
end
    ConstTerm = garchget(Coefficients, 'C');

if p > 0
    ARcoefSEE = Errors.AR;
end

if q > 0
    MAcoefSEE = Errors.MA;
end
ConstTermSEE = Errors.C;

%checking convergance
IsConverged = strcmp(Summery.converge, 'Function converged to a solution');
%checking stationarity of parameters
IsStatInvert = strcmp(Summery.warning, 'No warnings'); 
end


function [Spec]= SelectBestModel()
global df teta
bestp=0; 
bestq=0;
bestConstTerm = 0;
bestARcoef = [];
bestMAcoef = [];
BestRMSE = 1000000000;

for i=0:2
    for j=0:2
        if i==0 && j==0
            continue;
        end
        
        [ARcoef, MAcoef, ConstTerm, ARcoefSEE, MAcoefSEE, ConstTermSEE, Residuals, IsStatInvert, IsConverged] = Estimate(df.finalSeries, i, j);
        if ~IsStatInvert || ~IsConverged
            continue;
        end
        
        if checkSignificance(ARcoef, MAcoef, ARcoefSEE, MAcoefSEE)==0
            continue;
        end
        
        if CheckResiduals(Residuals)==1
            continue;
        end
        
        RMSE = computeRMSE(Residuals, i, j);
        strRMSE = sprintf('\tBest RMSE=%f, \tRMSE=%f', BestRMSE, RMSE);
        disp(strRMSE);
        if(RMSE < BestRMSE)
            BestRMSE = RMSE;
            field1 = 'P';  
            value1 = i;
            field2 = 'Q';  
            value2 = j;
            field3 = 'ARcoef';  
            value3 = ARcoef;
            field4 = 'MAcoef';  
            value4 = MAcoef;
            field5 = 'ConstTerm';  
            value5 = ConstTerm;
            
            Spec = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5);
        end
    end
end
end


function [RMSE] = computeRMSE(Residuals, p, q)
SSE =0;
n= length(Residuals);
for i=1:n
    SSE = SSE + (Residuals(i)*Residuals(i));
end
RMSE = sqrt(SSE/(n-p-q));
end


function [result] = checkSignificance(ARcoef, MAcoef, ARcoefSEE, MAcoefSEE)
result = 1;
end


function [H] = CheckResiduals(Residuals)
    [H, pValue, W] = ShapiroWilkTest(Residuals, 0.0006, 1);
end
