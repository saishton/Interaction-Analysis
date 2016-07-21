function [mins,maxs] = mm_ExpMLGPWei(FT,MLE,Moments)

brokenMins = Inf(1,8);
brokenMaxs = -Inf(1,8);

FT_Vec = zeros(1,8);
MLE_Vec = zeros(1,8);
Moments_Vec = zeros(1,8);

FT_Vec(1) = FT.Exponential.Parameters.Scale;
FT_Vec(2) = FT.MittagLeffler.Parameters.Stability;
FT_Vec(3) = FT.MittagLeffler.Parameters.Scale;
FT_Vec(4) = FT.GenPareto.Parameters.Shape;
FT_Vec(5) = FT.GenPareto.Parameters.Scale;
FT_Vec(6) = FT.GenPareto.Parameters.Location;
FT_Vec(7) = FT.Weibull.Parameters.Scale;
FT_Vec(8) = FT.Weibull.Parameters.Shape;

MLE_Vec(1) = MLE.Exponential.Parameters.Scale;
MLE_Vec(2) = NaN;
MLE_Vec(3) = NaN;
MLE_Vec(4) = NaN;
MLE_Vec(5) = NaN;
MLE_Vec(6) = NaN;
MLE_Vec(7) = MLE.Weibull.Parameters.Scale;
MLE_Vec(8) = MLE.Weibull.Parameters.Shape;

Moments_Vec(1) = Moments.Exponential.Parameters.Scale;
Moments_Vec(2) = NaN;
Moments_Vec(3) = NaN;
Moments_Vec(4) = Moments.GenPareto.Parameters.Shape;
Moments_Vec(5) = Moments.GenPareto.Parameters.Scale;
Moments_Vec(6) = Moments.GenPareto.Parameters.Location;
Moments_Vec(7) = NaN;
Moments_Vec(8) = NaN;

mins = min([brokenMins;FT_Vec;MLE_Vec;Moments_Vec]);
maxs = max([brokenMaxs;FT_Vec;MLE_Vec;Moments_Vec]);