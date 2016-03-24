function fit = testStatistics(X,Z)

n = length(X);
one = ones(n,1);
iVec = 1:n;
iVec = iVec';

Dplu = max((iVec/n)-Z);
Dmin = max(Z-((iVec-one)/n));
D = max(Dplu,Dmin);

CvM_vec = (Z-(2*iVec-one)/(2*n)).^2;
Wsq = sum(CvM_vec)+(1/(12*n));

V = Dplu+Dmin;

WatMod = n*(sum(Z)/n - 0.5)^2;
Usq = Wsq-WatMod;

Zswitch = flipud(Z);
AD_vec = (2*iVec-one).*(log(Z)+log(1-Zswitch));
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