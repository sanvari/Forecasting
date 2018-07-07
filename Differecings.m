classdef Differecings < handle
    % All methods ralated to differencing used in the program are defined
    % here.
    
    properties (SetAccess = public)
        oldSeries;
        logSeries;
        normSeries;
        OriginalSeries;
        finalSeries;
        nDiff;
        lag =18;
        SeriesSize;
        NonSeasonalSeries;
        FinalForcastedSeries;
        FinalTransformedSeries;
        VarStatSeries;
        FinalSeasonalSeries;
        SeSeries;
        Sx;
        Sy;
    end
    
    methods
        function dff = Differecings(Series)
            dff.oldSeries = Series;
            [dff.SeriesSize, cols] = size(Series);
            dff.nDiff = 0;
        end
        
        
        function LogFiff(dff)
            %            dff.normSeries = dff.oldSeries;
            dff.logSeries = dff.NonSeasonalSeries;
            ss = length(dff.logSeries);
            flag = 1;
            
            [h1, pValue1, stat1] = kpsstest(dff.logSeries);
            
            for e=1:floor(ss/2)
                dff.Sx(e) = dff.logSeries(e);
            end
            for e=floor(ss/2)+1:ss
                dff.Sy(e)= dff.logSeries(e);
            end
            
            [ht2]= myTTest2(dff.Sx,dff.Sy);
            %    [ht,stat3,sizet] = myTTest(dff.logSeries, ss, flag);
            
            %    [h,p,Qstat,crit] = lbqtest(dff.logSeries,'Lags',[5,10,20])
            
            
            while h1 ==1 || ht2==1  || any(h)  %any(h)1==1
                [ss,cs] = size(dff.logSeries);
                dff.logSeries = diff(log(dff.logSeries));
                dff.nDiff = dff.nDiff + 1;
                
                for e=1:floor(ss/2)
                    dff.Sx(e) = dff.logSeries(e);
                end
                for e=floor(ss/2)+1:ss
                    dff.Sy(e)= dff.logSeries(e);
                end
                
                [ht2]= myTTest2(dff.Sx,dff.Sy);

                %                    [h1, pValue2, stat2] = kpsstest(dff.logSeries);
                [ht,stat3,sizet] = myTTest(dff.logSeries,ss,flag);
                [h,p,Qstat,crit] = lbqtest(dff.logSeries,'Lags',[5,10,20])
            end
            dff.finalSeries = dff.logSeries;
        end
        
        function LogTrans(dff)
            l = length(dff.oldSeries);
            for q=1:l
                dff.VarStatSeries(q,1) = log(dff.oldSeries(q,1)); 
            end
            figure(23);
            plot(dff.VarStatSeries);
            xlim([0, l]);
            set(gca,'XTick', 0:10:l);
            set(gca,'XTickLabel', 0:10:l);
            title('log Transformed Data');
        end
            
        function NormalDiff(dff)
            % lbqtest :  assesses the null hypothesis that a series of residuals exhibits no autocorrelation for a fixed
            %            number of lags L, against the alternative that some autocorrelation coefficient rho(k), k = 1, ..., L, is nonzero
            
             %The KPSS test assesses the null hypothesis that a univariate time series y is trend stationary
             %           against the alternative that it is a nonstationary unit-root process
                
          %  dff.normSeries = dff.oldSeries;
            dff.normSeries = dff.NonSeasonalSeries;
          %  dff.normSeries = dff.VarStatSeries;
            ss = length(dff.normSeries);
            flag = 1;
            
            half = floor(ss/2);
            dff.Sx = dff.normSeries(1:half);
            dff.Sy = dff.normSeries(half+1:ss);
            
            [ht2]= myTTest2(dff.Sx,dff.Sy);
            [h1, pValue1, stat1] = kpsstest(dff.normSeries);
            [ha,pValue,stat,cValue,reg] = adftest(dff.normSeries, 'alpha', 0.1);
            
   %         [h,p,Qstat,crit] = lbqtest(dff.normSeries, 'alpha', 0.001);
             [hc] = CorrelationTest(dff.normSeries);
             [acf,lags,bounds] = autocorr(dff.normSeries);
         %   [hr] = runstest(acf,mean(acf), 'alpha', 0.1);
          %   [h,p,Qstat,crit] = lbqtest(acf, 'alpha', 0.1);
          %  [hr] = runstest(dff.normSeries,mean(dff.normSeries), 'alpha', 0.1);
             [h,p,Qstat,crit] = lbqtest(acf,'Lags',[5 , 10, 15], 'alpha', 0.01);
            
            
            % if h1=0 --> OK (trend stationary)
            % if ht2=0 -->OK (mean stationary)
            % if ha=1 --> OK (is not unit root)
            % if hc=0 -->OK (is not correlated)
            % if hr=0 -->OK (sequence is random)
            
            while h1==1 || ht2==1  ||ha ==0 || hc==1 || any(h)==1 %|| hr==1%   
                dff.normSeries = diff(dff.normSeries);
                ss = length(dff.normSeries);
                dff.nDiff = dff.nDiff + 1;
                half = floor(ss/2);
                dff.Sx = dff.normSeries(1:half);
                dff.Sy = dff.normSeries(half+1:ss);
                
                [ht2]= myTTest2(dff.Sx,dff.Sy);
                [h1, pValue2, stat2] = kpsstest(dff.normSeries);
                [hc] = CorrelationTest(dff.normSeries);
                [ha,pValue,stat,cValue,reg] = adftest(dff.normSeries, 'alpha', 0.1);
                [acf,lags,bounds] = autocorr(dff.normSeries);
               % [hr] = runstest(dff.normSeries,mean(dff.normSeries), 'alpha', 0.001);
          
             %  [h,p,Qstat,crit] = lbqtest(acf, 'alpha', 0.01);
                [h,p,Qstat,crit] = lbqtest(acf,'Lags',[5 , 10, 15], 'alpha', 0.01);
         
             %   [hr] = runstest(acf,mean(acf), 'alpha', 0.001);
          
                
            end
            disp(dff.nDiff);
            dff.finalSeries = dff.normSeries;
            figure(13);
            plot(dff.finalSeries);
            xlim([0, ss]);
            set(gca,'XTick', 0:10:ss);
            set(gca,'XTickLabel', 0:10:ss);
            title('Stationary Data');
            
        end
        
        function seasonaldiff(dff)
            %(Series,lag)
            % This function performs seasonal differencing transformation on
            % input matrix.
            % Input:
            % Series: input matrix with m row and 1 column
            % Lag: seasonal differencing lag, for example 4,6,8,12,24
            % Output:
            % NonSeasonalSeries: output matrix withm?p rowand 1 column
            
            %df.SeSeries = df.oldSeries;
            dff.SeSeries = dff.VarStatSeries;
            [r,c]=size(dff.SeSeries);
            if dff.lag>0 && dff.lag < r
                for i=1:r-dff.lag
                    dff.NonSeasonalSeries(i,c)=dff.SeSeries(i+dff.lag)-dff.SeSeries(i);
                end
            else
                error('Invalid Lag!');
                dff.NonSeasonalSeries = dff.SeSeries;
            end
            figure(12);
            plot(dff.NonSeasonalSeries);
            ss2 = length(dff.NonSeasonalSeries);
            xlim([0, ss2]);
            set(gca,'XTick', 0:10:ss2);
            set(gca,'XTickLabel', 0:10:ss2);
            title('Nonseasonal Data');

        end
        
        function [NewSeries1]= inverseNormalDiff(dff, ForecastedSeries)
            % First determine how many numbers are there
            
            %[rows1, cols1]=size(dff.oldSeries);
            [rows1, cols1]=size(dff.NonSeasonalSeries);
            
            [rows2, cols2]=size(ForecastedSeries);
            
            % Put the original data
            %  NewSeries1 = dff.oldSeries;
            NewSeries1 = dff.NonSeasonalSeries;
            
            switch dff.nDiff
                case 0
                    for i=1:rows2
                        NewSeries1(i+rows1, cols1) = ForecastedSeries(i, cols2);
                    end
                    
                case 1
                    for i=1:rows2
                        NewSeries1(i+rows1, cols1) = NewSeries1(i+rows1-1, cols1) + ForecastedSeries(i, cols2);
                    end
                case 2
                    for i=1:rows2
                        NewSeries1(i+rows1, cols1) = ForecastedSeries(i, cols2) + 2*NewSeries1(i+rows1-1, cols1) - NewSeries1(i+rows1-2, cols1);
                    end
                case 3
                    for i=1:rows2
                        NewSeries1(i+rows1,cols1)= ForecastedSeries(i,cols2)+3*NewSeries1(i+rows1-1,cols1)-3*NewSeries1(i+rows1-2,cols1)+NewSeries1(i+rows1-3,cols1);
                    end
                otherwise
                    error('it is not good to have to process series that need more than 3 times differntiating');
            end

        end
        
        function [OriginalSeries]= inverseLogDiff(dff, ForecastedSeries)
            %        [r1, c1]=size(dff.oldSeries);
            [r1, c1]=size(dff.NonSeasonalSeries);
            
            [r2, c2]=size(ForecastedSeries);
            
            % Put the original data
            %      OriginalSeries = dff.oldSeries;
            OriginalSeries = dff.NonSeasonalSeries;
            
            switch dff.nDiff
                case 0
                    for i=1:r2
                        OriginalSeries(i+r1, c1) = ForecastedSeries(i, c2);
                    end
                case 1
                    for i=1:r2
                        OriginalSeries(i+r1, c1) = exp(log(OriginalSeries(i+r1-1, c1)) + ForecastedSeries(i, c2));
                    end
                    %     case 2
                    %        for i=1:rows2
                    %           NewSeries1(i+rows1, cols1) = ForecastedSeries(i, cols2) + 2*NewSeries1(i+rows1-1, cols1) - NewSeries1(i+rows1-2, cols1);
                    %      end
                    %   case 3
                    %      for i=1:rows2
                    %         NewSeries1(i+rows1,cols1)= ForecastedSeries(i,cols2)+3*NewSeries1(i+rows1-1,cols1)-3*NewSeries1(i+rows1-2,cols1)+NewSeries1(i+rows1-3,cols1);
                    %    end
                otherwise
                    error('it is not good to have to process series that need more than 3 times differntiating');
            end
            
        end
        
        function [FinalForcastedSeries]=invseasonaldiff(oldSeries,OriginalSeries,Lag)
            % This functionperforms reverse seasonal differencing on input series
            % Input:
            % SeriesBeforeSeasonalDiff: series before seasonal differencing
            % NonSeasonalSeries: seasonal differenced series
            % Lag: seasonal differencing lag
            % Output:
            % NewSeries:reverse seasonal differenced series
            [r1,c1]=size(oldSeries);
            [r2,c2]=size(OriginalSeries);
            if Lag==0
                FinalForcastedSeries=OriginalSeries;
            else
                point=r1-Lag;
                FinalForcastedSeries=oldSeries;
                for i=1:r2-point
                    FinalForcastedSeries(r1+i,c1)=FinalForcastedSeries(r1+i-Lag,c1)+OriginalSeries(point+i,c2);
                end
            end
            
            
        end
        
        function [FinalSeasonalSeries]=invseasonaldiffN(dff, NewSeries1)
            % This functionperforms reverse seasonal differencing on input series
            % Input:
            % SeriesBeforeSeasonalDiff: series before seasonal differencing
            % NonSeasonalSeries: seasonal differenced series
            % Lag: seasonal differencing lag
            % Output:
            % NewSeries:reverse seasonal differenced series
            
            %[r1,c1]=size(dff.oldSeries);
            [r1,c1]=size(dff.VarStatSeries);
            [r2,c2]=size(NewSeries1);
            Lag=dff.lag
            if Lag==0
                FinalSeasonalSeries=NewSeries1;
            else
                point=r1-Lag;
                %FinalSeasonalSeries=dff.oldSeries;
                FinalSeasonalSeries=dff.VarStatSeries;
                for i=1:r2-point
                    FinalSeasonalSeries(r1+i,c1)=FinalSeasonalSeries(r1+i-Lag,c1)+NewSeries1(point+i,c2);
                end
            end
        end
   
        function [FinalTransformedSeries] = inverseLogTransform(dff, FinalSeasonalSeries)
            [r4,c4]=size(FinalSeasonalSeries);
            for i=1:r4
                FinalTransformedSeries(i,1) = exp(FinalSeasonalSeries(i,1));
            end
        end
    end
    
end

