% Add path
addpath('FunctionFolder');
% definition of the basic parameters
ns=0.158453 + 3.56052*1i;
nb=1;
lamda=655e-9;
RAu=40e-9;
distance=RAu+14e-9;
nP=Fun_nP_DipoleNearSphere(lamda,nb,ns,RAu,distance);