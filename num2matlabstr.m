function [str] = num2matlabstr(r)

if r==Inf
    str = '\infty';
elseif r==-Inf
    str = '-\infty';
elseif isnan(r)
    str = 'X';
elseif r==0
    str = '0';
else
    expon = floor(log10(abs(r)));
    coeff = r/(10^expon);
    exponstr = num2str(expon);
    coeffstr = num2str(coeff,'%.4f');
    str = strcat(coeffstr,'\times10^{',exponstr,'}');
end
end