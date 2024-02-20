function [diwout_dt] = exoBootsODEfun(t, iwout, varargin)
    %Constants
    if isempty(varargin)
        R=1.03; %Units: Ohms
    else
        R=varargin{1};
    end 
    L=2.04e-4; %Units: H
    Kt=4.48e-2; %Units: Nm/A
    Kv=7.1*pi; %Units: rad/Vs
    Jm=1.01e-5; %Units: kg/m^2
    N=1; %Unitless
    bm=0.001; %Units: Nms/rad

    %Creating matrix and solving for derivatives
    a11=-R/L;
    a12=-N/Kv/L;
    a21=Kt/Jm/N;
    a22=-bm/Jm;
    A = [a11, a12; a21, a22];
    B = iwout;
    diwout_dt = A*B;
end