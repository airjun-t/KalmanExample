function [x_guess,  P_guess] = project(F,B, x,u,P, Q)
% best guess before filtering and considering the measurement
    x_guess = F * x + B * u;

    P_guess = F * P * F' + Q;
     
end