function fit = testStatistics(X,Z)

n = length(X);
uni = unique(X);
if length(uni)==n
    fullecdf = 0:1/n:1;
else
    h = hist(X,X);
    cumh = cumsum(h);
    fullecdf = [0 cumh/n];
end
fullecdf = fullecdf';
lowerecdf = fullecdf(1:n);
upperecdf = fullecdf(2:n+1);
middlecdf = (lowerecdf+upperecdf)/2;

Dplu = max(lowerecdf-Z);
Dmin = max(Z-upperecdf);
D = max(Dplu,Dmin);

CvM_vec = (Z-middlecdf).^2;
Wsq = sum(CvM_vec)+(1/(12*n));

V = Dplu+Dmin;

WatMod = n*(mean(Z) - 0.5)^2;
Usq = abs(Wsq-WatMod);

Zswitch = flipud(Z);
AD_vec = (2*middlecdf).*(log(Z)+log(1-Zswitch));
Asq = -sum(AD_vec) - n;
if Asq == inf
    Asq = 10^50;
elseif Asq == -inf
    Asq = -10^50;
end

fit = struct(   'Kolmogorov_D',D,...
                'Cramer_von_Mises',Wsq,...
                'Kuiper',V,...
                'Watson',Usq,...
                'Anderson_Darling',Asq);
end