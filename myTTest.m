            %%%%%%%%%%%%%%%%%%%%%%%%%%%   t-test Begin     %%%%%%%%%%%%%%%%%%55
        function [ht,stats,sizeT]=myTTest(x,n,flag)
            % This function calculates t-test statistics for input matrix x and
            % performs t-test hypothesis test.
            % Input:
            % x: input matrix
            % n: total number of original series data
            % flag: 0 for ACF and 1 for PACF
            % Output:
            % h: binary value for the result of t-test. 0 pass, 1 not pass.
            % stats: t-test statistics
            % size: total number of t-test statistics
            if nargin <1,
                error('Requires at least one input argument.');
            end
            [m1 , n1]=size(x);
            if (m1 ~=1 && n1 ~=1)
                error('First argument has to be a vector.');
            end
            x=x(~isnan(x));
            sizeT=length(x);
            if flag==0
                rootsum=zeros(sizeT,1);
                for i=2:sizeT
                    for j=1:i-1
                        rootsum(i)=rootsum(i)+x(j)^2;
                    end
                end
                stats=x.*sqrt(n)./sqrt(1+2.*rootsum);
            else
                stats=x.*sqrt(n);
            end
            % Rule of Thumb for critical level of t-test
            crit=2;
            % Determine if the actual significance exceeds the desired significance
            ht=0;
            n=0;
            for i=2:sizeT
                if abs(stats(i))>=crit
                    n=n+1;
                end
            end
            if n>0
                ht=1;
            end
        end

