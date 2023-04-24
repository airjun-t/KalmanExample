function [new_state,state_dot] = rk4(func_handle, time, dt, state, inputs)
%RK4 Runge-Kutta 4 Matlab implementation
%   Detailed explanation goes here
    k1 = func_handle(time     , state            , inputs);
    k2 = func_handle(time+dt/2, state + dt*0.5*k1, inputs);
    k3 = func_handle(time+dt/2, state + dt*0.5*k2, inputs);
    k4 = func_handle(time+dt  , state + dt*k3    , inputs);

    state_dot = (k1 + 2.*(k2 + k3) + k4)./6.0;
    new_state = state + dt.*state_dot;
end
