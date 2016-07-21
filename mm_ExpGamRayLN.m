function [mins,maxs] = mm_ExpGamRayLN(FT,MLE,Moments)

brokenMins = Inf(1,6);
brokenMaxs = -Inf(1,6);

FT_Vec = zeros(1,6);
MLE_Vec = zeros(1,6);
Moments_Vec = zeros(1,6);

FT_Vec(1) = FT.Exponential.Parameters.Scale;
FT_Vec(2) = FT.Gamma.Parameters.Shape;
FT_Vec(3) = FT.Gamma.Parameters.Scale;
FT_Vec(4) = FT.Rayleigh.Parameters.Scale;
FT_Vec(5) = FT.LogNormal.Parameters.Location;
FT_Vec(6) = FT.LogNormal.Parameters.Scale;

MLE_Vec(1) = MLE.Exponential.Parameters.Scale;
MLE_Vec(2) = MLE.Gamma.Parameters.Shape;
MLE_Vec(3) = MLE.Gamma.Parameters.Scale;
MLE_Vec(4) = MLE.Rayleigh.Parameters.Scale;
MLE_Vec(5) = MLE.LogNormal.Parameters.Location;
MLE_Vec(6) = MLE.LogNormal.Parameters.Scale;

Moments_Vec(1) = Moments.Exponential.Parameters.Scale;
Moments_Vec(2) = Moments.Gamma.Parameters.Shape;
Moments_Vec(3) = Moments.Gamma.Parameters.Scale;
Moments_Vec(4) = Moments.Rayleigh.Parameters.Scale;
Moments_Vec(5) = Moments.LogNormal.Parameters.Location;
Moments_Vec(6) = Moments.LogNormal.Parameters.Scale;


mins = min([brokenMins;FT_Vec;MLE_Vec;Moments_Vec]);
maxs = max([brokenMaxs;FT_Vec;MLE_Vec;Moments_Vec]);