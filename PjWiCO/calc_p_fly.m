function P = calc_p_fly(v)
%   CALC_P_PFLY calculates Power based on power consumption model
%   P = calc_p_fly(v) returns power depending on velocity of UAV v
    P_o = 80; P_i = 88.6; Sigma = 300; r=0.4; v_o = 4.03; % UAV parameters of power consumtption model
    P = P_o * (3*norm(v)^2 / (Sigma^2*r^2)) ...
        + P_i * (sqrt(sqrt(1+ norm(v)^4/ (4*v_o^4))- norm(v)^2 / (2*v_o^2))-1);
end