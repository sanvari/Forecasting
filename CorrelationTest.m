 function [hc] = CorrelationTest(y)
 
 [acf,lags,bounds] = autocorr(y);
% [pacf,lags,bounds] = parcorr(y);
 sy = length(acf);
 %sy = length(pacf);
 
 hc = 0;
 sig = 0;
 for i=1:10
     x = acf(i, 1);
   %  if x < bounds(2,1) || x > bounds(1,1)
     if x < -0.2 || x > 0.2
         sig = sig+1;
     end
 end 

 if(sig>5)
     hc = 1;
end