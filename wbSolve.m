function [a,b] = wbSolve(M1,M2)

fnK = (M1^2)/M2;

syms k
eqn = (gamma(1+(1/k))^2)/(gamma(1+(2/k))) == fnK;
numK = vpasolve(eqn,k);
numLam = M1/(gamma(1+(1/numK)));

a=double(numLam);
b=double(numK);
end