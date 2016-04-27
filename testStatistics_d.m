function fit = testStatistics_d(X,Z,step)

n = length(X);
m = length(Z);
expected = zeros(1,m-1);

parfor i=1:m-1
    expected(i) = Z(i+1)-Z(i);
end
expected = expected*n;

grid = (min(X)-step):step:max(X);
observed = zeros(1,m-1);
parfor i=1:m-1
    observed(i) = numel(X( X(:)>grid(i) & X(:)<=grid(i+1)));
end

error = observed-expected;
errorsum = cumsum(error);
S = max(abs(errorsum));

chiSq = sum((error.^2)./expected);

fit = struct(   'Modified_Kol',S,...
                'Chi_Squared',chiSq,...
                'Root_MSE',0,...
                'R_Squared',0);
end