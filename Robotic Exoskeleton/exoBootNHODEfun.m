function [diwtheta_dt] = exoBootNHODEfun(t, iwtheta, N, exoTorque, exoAngle, exoTime)
%Constants
    R=1.03; %Units: Ohms
    L=2.04e-4; %Units: H
    Kt=4.48e-2; %Units: Nm/A
    Kv=7.1*pi; %Units: rad/Vs
    Jm=1.01e-5; %Units: kg/m^2
    bm=0.001; %Units: Nms/rad

    %Interpolating our data
    T_out = interp1(exoTime,exoTorque,t);
    theta_des = interp1(exoTime,exoAngle,t);
    theta_out = iwtheta(3,1);

    %Setting the value of e_in
    e_in = 5*(theta_des-theta_out)*N;
    if e_in > 48
        e_in = 48;
    elseif e_in < -48
        e_in = -48;
    end

    %Creating matrix and solving for derivatives
    a11=-R/L;
    a12=-N/Kv/L;
    a13=0;
    a21=Kt/Jm/N;
    a22=-bm/Jm;
    a23=0;
    a31=0;
    a32=1;
    a33=0;
    A = [a11, a12, a13; a21, a22, a23; a31, a32, a33];
    B = iwtheta;
    c1 = e_in/L;
    c2 = -T_out/Jm/N^2;
    c3 = 0;
    C=[c1; c2; c3];
    
    diwtheta_dt = A*B+C;
end