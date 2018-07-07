
%%%%%%%%%%%%%%%%%%%%     shapiro-wilk Begin    %%%%%%%%%%%%%

function [H, pValue, W] = ShapiroWilkTest(x, alpha, tail)

% Ensure the sample data is a VECTOR.
if numel(x) == length(x)
    x  =  x(:);               % Ensure a column vector.
else
    error(' Input sample ''X'' must be a vector.');
end

% Remove missing observations indicated by NaN's
x  =  x(~isnan(x));

%checking sample size.
if length(x) < 3
    error(' Sample vector ''X'' must have at least 3 valid observations.');
end
if length(x) > 5000
    warning('Shapiro-Wilk test might be inaccurate due to large sample size ( > 5000).');
end

if nargin < 1,
    error('Requires at least one input argument.');
    return;
end
% Ensure the significance level, ALPHA, is a % scalar, and set default if necessary.
if (nargin >= 2) && ~isempty(alpha)
    if numel(alpha) > 1
        error(' Significance level ''Alpha'' must be a scalar.');
        return;
    end
    if (alpha <= 0 || alpha >= 1)
        error(' Significance level ''Alpha'' must be between 0 and 1.');
        return;
    end
else
    alpha  =  0.05;
end

% Ensure the type-of-test indicator, TAIL, is a scalar integer from
% the allowable set {-1 , 0 , 1}, and set default if necessary.
if (nargin >= 3) && ~isempty(tail)
    if numel(tail) > 1
        error('Type-of-test indicator ''Tail'' must be a scalar.');
        return;
    end
    if (tail ~= -1) && (tail ~= 0) && (tail ~= 1)
        error('Type-of-test indicator ''Tail'' must be -1, 0, or 1.');
        return;
    end
else
    tail  =  1;
end
% calculating a's for weights as a function of the m's
% See Royston (1995) for details in the approximation.
x       =   sort(x); % Sort the vector X in ascending order.
n       =   length(x);
m       =   norminv(((1:n)' - 3/8) / (n + 0.25));
weights =   zeros(n,1); % Preallocate the weights.


if kurtosis(x) > 3
    % The Shapiro-Francia test is better for leptokurtic samples.
    weights =   1/sqrt(m'*m) * m;
    W   =   (weights' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
    
    % The Shapiro-Francia statistic W is calculated to avoid excessive rounding
    % errors for W close to 1 (a potential problem in very large samples).
    
    nu      =   log(n);
    u1      =   log(nu) - nu;
    u2      =   log(nu) + 2/nu;
    mu      =   -1.2725 + (1.0521 * u1);
    sigma   =   1.0308 - (0.26758 * u2);
    
    newSFstatistic  =   log(1 - W);
    
    % Compute the normalized Shapiro-Francia statistic and its p-value.
    NormalSFstatistic =   (newSFstatistic - mu) / sigma;
    
    % the next p-value is for the tail = 1 test.
    pValue   =   1 - normcdf(NormalSFstatistic, 0, 1);
else
    % The Shapiro-Wilk test is better for platykurtic samples.
    c    =   1/sqrt(m'*m) * m;
    u    =   1/sqrt(n);
    p1   =   [-2.706056,4.434685,-2.071190,-0.147981,0.221157,c(n)];
    p2   =   [-3.582633,5.682633,-1.752461,-0.293762,0.042981,c(n-1)];
    
    p3   =   [-0.0006714 , 0.0250540 , -0.39978 , 0.54400];
    p4   =   [-0.0020322 , 0.0627670 , -0.77857 , 1.38220];
    p5   =   [0.00389150 , -0.083751 , -0.31082 , -1.5861];
    p6   =   [0.00303020 , -0.082676 , -0.48030];
    
    p7   =   [0.459 , -2.273];
    
    weights(n)   =   polyval(p1 , u);
    weights(1)   =   -weights(n);
    
    % Special case when n=3
    if n == 3
        weights(1)  =   0.707106781;
        weights(n)  =   -weights(1);
    end
    
    if n >= 6
        weights(n-1) =   polyval(p2 , u);
        weights(2)   =   -weights(n-1);
        
        count  =   3;
        phi    =   (m'*m - 2 * m(n)^2 - 2 * m(n-1)^2) / ...
            (1 - 2 * weights(n)^2 - 2 * weights(n-1)^2);
    else
        count  =   2;
        phi    =   (m'*m - 2 * m(n)^2) /(1 - 2 * weights(n)^2);
    end
    
    % The vector 'WEIGHTS' obtained next corresponds to the same coefficients
    % listed by Shapiro-Wilk in their original test for small samples.
    weights(count : n-count+1)  =  m(count : n-count+1) / sqrt(phi);
    
    % The Shapiro-Wilk statistic W is calculated to avoid excessive rounding
    % errors for W close to 1 (a potential problem in very large samples).
    W   =   (weights' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
    
    % Calculate the significance level for W (exact for n=3).
    newn    =   log(n);
    
    if (n < 3)
        error('n is too small.');
        return;
    end
    if (n > 3) && (n <= 11)
        mu      =   polyval(p3 , n);
        sigma   =   exp(polyval(p4 , n));
        gam     =   polyval(p7 , n);
        newSWstatistic  =   -log(gam-log(1-W));
        
    elseif (n >= 12) && (n <=5000)
        
        mu      =   polyval(p5 , newn);
        sigma   =   exp(polyval(p6 , newn));
        gam=0;
        newSWstatistic  =   log(1 - W)+gam;
        
    elseif n == 3
        mu      =   0;
        sigma   =   1;
        newSWstatistic  =   0;
    else
        error('n is not in the proper size range.'); %error('n is too large.');return,
    end
    
    % Compute the normalized Shapiro-Wilk statistic and its p-value.
    NormalSWstatistic       =   (newSWstatistic - mu) / sigma;
    
    % The next p-value is for the tail = 1 test.
    pValue       =   1 - normcdf(NormalSWstatistic, 0, 1);
    
    % Special attention when n=3 (this is a special case).
    if n == 3
        pValue  =   1.909859 * (asin(sqrt(W)) - 1.047198);
        NormalSWstatistic =   norminv(pValue, 0, 1);
    end
end

% The p-value just found is for the tail = 1 test.
if tail == 0
    pValue = 2 * min(pValue, 1-pValue);
elseif tail == -1
    pValue = 1 - pValue;
end
% To maintain consistency with existing Statistics Toolbox hypothesis
% tests, returning 'H = 0' implies that we 'Do not reject the null
% hypothesis at the significance level of alpha' and 'H = 1' implies
% that we 'Reject the null hypothesis at significance level of alpha.'
% Null Hypothesis: X is normal with unspecified mean and variance
if pValue >= alpha;
    H=0;
    disp('Data analyzed have a normal distribution.');
else
    disp('Data analyzed do not have a normal distribution.');
    H=1;
end
end
%%%%%%%%%%%%%%%%%%%%     shapiro-wilk end    %%%%%%%%%%%%%

