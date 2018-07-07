            %%%%%%%%%%%%%%%%%%%%%%%%%%%   t-test Begin     %%%%%%%%%%%%%%%%%%55
        function [ht2]=myTTest2(x,y)
            % This function decide for  mean stationarit of series x and y and
            % performs t-test hypothesis test.
            % Input:
            % x: input matrix
            % y: input matrix
            % nx: total number of original series x
            % ny: total number of original series y
            
            % Output:
            % ht2: binary value for the result of t-test. 0 pass, 1 not pass.

            if nargin <1,
                error('Requires at least one input argument.');
            end
            [m1 , n1]=size(x);
            if (m1 ~=1 && n1 ~=1)
                error('First argument has to be a vector.');
            end
            [m2 , n2]=size(y);
            if (m2 ~=1 && n2 ~=1)
                error('First argument has to be a vector.');
            end
            x=x(~isnan(x));
            y=y(~isnan(y));
            sizex=length(x);
            sizey=length(y);
            Mx = mean(x);
            My = mean(y);
           
            t = (Mx - My)/sqrt((var(x)/sizex)+(var(y)/sizey));
            if t >= 2
                ht2 = 1;
            else
                ht2 = 0;
            end
      
        end

