classdef MyForecast
    %UNTÝTLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        EndedDif;
        z = 0;
    end
    
    methods
        function frc = MyForecast(df)
            frc.EndedDif = df;
        end
        
        function [RMSE] = computeRMSE(frc, Residuals, p, q)
            SSE =0;
            n= length(Residuals);
            for i=1:n
                SSE = SSE + (Residuals(i)*Residuals(i));
            end
            RMSE = sqrt(SSE/(n-p-q));
        end
        
        function [result] = checkSignificance(frc, ARcoef, MAcoef, ARcoefSEE, MAcoefSEE)
            result = 0;
            [s1,c1]= size(ARcoef);
            [s2,c2]= size(MAcoef);
            for i=1:s1
                if abs((ARcoef(i)/ARcoefSEE(i))) < 2
                    result = 1;
                    continue;
                end
            end
            for i=1:s2
                if abs((MAcoef(i)/MAcoefSEE(i))) < 2
                    result = 1;
                    continue;
                end
            end
        end
        
        function [H] = CheckNormalResiduals(frc, Residuals)
              [s2, c1] = size(Residuals);
         %   if s2<= 5000
          %        figure(15);
          %        qqplot(Residuals);
                  [H1, pValue, W] = ShapiroWilkTest(Residuals, 0.001, 1);
                % small alpha = easier accept of normal dist
       
        %          % hj = 0 data is normal
         %         % hj = 1 data is not normal
               %    hj = ttest(Residuals);
         
         %     else
         
           %     figure(15);
           %     qqplot(Residuals);
           
          %      normplot(Residuals);
                hj = jbtest(Residuals, 0.01);
           %    hj = lillietest(Residuals);
              % kstest
               
        %      end
                 if H1 == 0 || hj == 0
                      H = 0;
                 else
                      H =1;
                 end
        end
        
        %This function checks whether the residuals are random or there is
        %correlation between them.
        % test of Ljung and Box assesses the null hypothesis that a series of residuals exhibits no autocorrelation
        %for a fixed number of lags L, against the alternative that some autocorrelation coefficient rho(k),
        %k = 1, ..., L, is nonzero
        % If h=0 then no correlation is in residulas or they are random
        function [h, p, Qstat, crit]= checkRandomResidulas(frc, Residuals, i, j)
        %    [acf,lags,bounds] = autocorr(Residuals);
            [h,p,Qstat,crit] = lbqtest(Residuals,'Lags',[5 , 10, 15],'alpha',0.01) %,'dof',[5-i-j,10-i-j,15-i-j]);
        end
        
        
        function [ARcoef, MAcoef, ConstTerm, ARcoefSEE, MAcoefSEE, ConstTermSEE, Residuals, IsStatInvert, IsConverged, Sig] = Estimate(frc, p, q)
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
            ConstTerm = mean(frc.EndedDif.finalSeries)*(1-SumOfARCoef);
            
            %Model specification generation
            ModelSpec = garchset('C', ConstTerm, 'R', p, 'M',q, 'AR', ARcoef, 'MA', MAcoef, 'Display','off');
            
            %Coefficient Estimation
            % Coefficients - ARMAX/GARCH model specification structure containing the parameter estimates. Coefficients has the same structure as Spec.
 
            % Errors - ARMAX/GARCH model specification structure containing the standard errors of the parameter estimates. Errors has the same structure as Spec.
 
            %LLF - Optimized loglikelihood objective function value associated with the parameter estimates found in Coefficients.
 
            % Residuals - Vector of inferred innovations (fit residuals), the same size as finalSeries.
 
            % sigma - Vector of conditional standard deviations of Residuals, the same size as Residuals.
 
            %summary - Structure of summary information for the maximum likelihood estimation, with the following fields and possible values.
            [Coefficients, Errors, LLF,Residuals, Sigma, Summery]= garchfit(ModelSpec, frc.EndedDif.finalSeries);
            %            disp(Coefficients);
            %            disp(Errors);
            %            disp('Residuals');
            %            disp(Residuals);
            
            %extracting values
            Sig = Sigma;
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
        
        function [Spec]= SelectBestModel(frc)
            bestp=0;
            bestq=0;
            bestConstTerm = 0;
            bestARcoef = [];
            bestMAcoef = [];
            BestRMSE = 1000000000;
            
            for i=0:3
                for j=0:3
                    if i==0 && j==0
                        continue;
                    end
                    
                    [ARcoef, MAcoef, ConstTerm, ARcoefSEE, MAcoefSEE, ConstTermSEE, Residuals, IsStatInvert, IsConverged, Sig] = frc.Estimate(i, j);
                    if ~IsStatInvert || ~IsConverged
                        continue;
                    end
                    % checking significancy : all estimated parameters should be significantly different
                    %than zero(t-ratio significant)
                    if frc.checkSignificance(ARcoef, MAcoef, ARcoefSEE, MAcoefSEE) == 1
                        continue;
                    end
                    %checking stationarity of AR
                    stat = 0;
                    [s1,c1] =  size(ARcoef);
                    for n=1:s1
                        stat =  stat + ARcoef(n);
                    end
                    if stat>1
                        continue;
                    end
                    
                    %checking Invertibility of MA
                    Inv = 0;
                    [I1,c1] =  size(MAcoef);
                    for m=1:I1
                        Inv =  Inv + MAcoef(m);
                    end
                    if Inv>1
                        continue;
                    end
             
                      %checking residuals in ACF and PACF
                    figure(24);
                    subplot(2,1,1);
                    autocorr(Residuals);
                    subplot(2,1,2);
                    parcorr(Residuals);
                    title('ACF and PACF of Residuals');
                    p=0;
                    [sr,cr] =  size(Residuals);
                    level = 2/sqrt(sr);
                    
                    [acf,lags,bounds]= autocorr(Residuals);
                    [sf,cr] =  size(acf);
                    for f=1:sf
                        if acf(f)> level || acf(f)<-level
                            p = p+1;
                        end
                    end
                    if p > sf/4      %ceil(sf/4)
                        continue;
                    end
                    
                    [h,p,Qstat,crit] = frc.checkRandomResidulas(Residuals, i, j);
                    NoRand = 0;
                    nh = numel(h);
                    for a=1:nh
                        if h(a)==1
                            NoRand=NoRand+1
                        end
                    end
                    if  NoRand> 1
               %    if  any(h) ==1
                        disp('Residuals are not random!');
                  %      continue;
                    end
                    
                    % % checking randomness of residuals % %%%%%%%%%%%%%%
                    % in runstest, if h=0 , residuals are random
                   [h] = runstest(Residuals);
                   if h==1
                       continue;
                   end
                    
                   %checking heteroscedasticity of Residuals
       %            h = archtest(Residuals)
    %               if h==1
    %                  continue;
     %              end
                   % checking Normality of Residuals
                    if frc.CheckNormalResiduals(Residuals)==1
                        continue;
                    end
                    
                    
                    RMSE = frc.computeRMSE(Residuals, i, j);
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
                        field6 = 'Residuals';
                        value6 = Residuals;
                        
                        Spec = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5, field6, value6);
                    end
                end
            end
        end
        
        function DoForecast(frc,df)
            [Spec] = frc.SelectBestModel();
            ROrder = Spec.P;
            MOrder = Spec.Q;
            I = df.nDiff +1;
            ARcoef = Spec.ARcoef;
            MAcoef = Spec.MAcoef;
            ConstTerm = Spec.ConstTerm;
            % disp(frc.EndedDif.finalSeries);
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
            
            NumPeriods = 5;
            %%%%%%%%%%%%%%%%%%%%New%%%%%%%NEW%%%%%%%%%%%%%Version%%%%%%%NEW
            
            % model = arima(ROrder, frc.EndedDif.nDiff, MOrder);
            % model.Constant = ConstTerm; %NaN
            % for i = 1:ROrder
            %     model.AR{i:i} = ARcoef(1, i);
            %end
            
            %for i = 1:MOrder
            %  model.MA{i:i} = MAcoef(1, i);
            %end
            
            
            %model.Variance = sqrt(var(Spec.Residuals));
            
            %%%%model = arima('Constant', ConstTerm, 'AR', ARcoef, 'MA', MAcoef, 'Variance', V);
            %%%%model = arima('Constant', ConstTerm, 'AR', ARcoef, 'MA', MAcoef, 'Variance', V, 'P', ROrder, 'D', frc.EndedDif.nDiff, 'Q', MOrder);
            %%%%%fit = estimate(model, Y);
            
            %[X, YMSE] = forecast(model, NumPeriods);
            %          %%%%% fit = estimate(model, frc.EndedDif.finalSeries);
            %         %%%%% [X2, YMSE2] = forecast(fit, NumPeriods, 'Y0', frc.EndedDif.finalSeries);
            % newData = [frc.EndedDif.finalSeries; YMSE];
            % figure(10);
            % plot(newData );
            
            % disp('MeanForecast');
            % disp(YMSE);
            
            % disp('MeanRMSE');
            % disp(X);
            
            
            %%%%%%%%%%%%%%%%%%%%New%%%%%%%%%%%%%%%%%%%%Version%%%%%%%New
            
            %%%%%%%%%%%%%%%%%%%%Old Version%%%%%%%Old
            %old estimation
            NewSpec = garchset('R',ROrder,'M',MOrder,'C',ConstTerm,'AR',ARcoef, 'MA',MAcoef,'K', 0.001,'P', 0,'Q',0,'Display', 'On');
            [SigmaForecast, MeanForecast, SigmaTotal, MeanRMSE]= garchpred(NewSpec, frc.EndedDif.finalSeries, NumPeriods);
            
            %%[EstSpec, EstSE] = garchfit(NewSpec, frc.EndedDif.finalSeries);
            %%[SigmaForecast, MeanForecast, SigmaTotal, MeanRMSE]= garchpred(EstSpec, frc.EndedDif.finalSeries, NumPeriods);
            
            Newdata = [frc.EndedDif.finalSeries; MeanForecast];
            
            
            figure(10);
            plot(Newdata);
            title('DifferencedForecastedData');
            
            disp('MeanForecast');
            disp(MeanForecast);
            
            disp('MeanRMSE');
            disp(MeanRMSE);
            
            disp('SigmaForecast');
            disp(SigmaForecast);
            
            disp('SigmaTotal');
            disp(SigmaTotal);
            
            %%%%%%%%%%%%%%%%%%%%Old Version%%%%%%%Old
            
            %      [newSeries]= frc.EndedDif.inverseNormalDiff(YMSE);
            %            [newSeries]= frc.EndedDif.inverseLogDiff(YMSE);
            
            [newSeries]= frc.EndedDif.inverseNormalDiff(MeanForecast);
            %    [newSeries]= frc.EndedDif.inverseLogDiff(MeanForecast);
            %[Finalforcasted]=frc.EndedDif.invseasonaldiff(newSeries);
            
            [Sesonalforcasted]=frc.EndedDif.invseasonaldiffN(newSeries);
            
            [Finalforcasted]=frc.EndedDif.inverseLogTransform(Sesonalforcasted);
                      
            [s1, c1] = size(Finalforcasted);
         
            disp('period forecasted');
            for i=0:NumPeriods
                disp(Finalforcasted(s1-NumPeriods+i, c1));
            end
            figure(11)
            plot(Finalforcasted)
            xlim([0, s1])
            set(gca,'XTick', 0:10:s1)
            set(gca,'XTickLabel', 0:10:s1)
            title('EndedForecastedData');
        end
    end
end