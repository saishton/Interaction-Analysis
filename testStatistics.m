function fit = testStatistics(X,Z)

n = length(X);
fullecdf = zeros(n+1,1);
parfor i=1:n
   thisX = X(i);
   fullecdf(i+1) =  sum(X <= thisX);
end
lowerecdf = fullecdf(1:n);
upperecdf = fullecdf(2:n+1);
middlecdf = (lowerecdf+upperecdf)/2;

Dplu = max((lowerecdf/n)-Z);
Dmin = max(Z-(upperecdf/n));
D = max(Dplu,Dmin);

CvM_vec = (Z-middlecdf).^2;
Wsq = sum(CvM_vec)+(1/(12*n));

V = Dplu+Dmin;

WatMod = n*(sum(Z)/n - 0.5)^2;
Usq = Wsq-WatMod;

Zswitch = flipud(Z);
AD_vec = (2*middlecdf).*(log(Z)+log(1-Zswitch));
Asq = -sum(AD_vec)/n - n;

fit = struct(   'Kolmogorov_D_Plus',Dplu,...
                'Kolmogorov_D_Minus',Dmin,...
                'Kolmogorov_D',D,...
                'Cramer_von_Mises',Wsq,...
                'Kuiper',V,...
                'Watson',Usq,...
                'Anderson_Darling',Asq,...
                'Root_MSE',0,...
                'R_Squared',0);
end