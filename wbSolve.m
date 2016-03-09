function [a,b] = wbSolve(M1,M2)

syms lambda k
eqn1 = lambda*gamma(1+(1/k)) == M1;
eqn2 = lambda^2*gamma(1+(2/k)) == M2;
eqn = [eqn1, eqn2];
para = [lambda, k];
[lambdaSol,kSol] = vpasolve(eqn,para);
 
a=double(lambdaSol);
b=double(kSol);
end