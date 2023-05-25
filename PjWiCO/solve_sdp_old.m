function optim_var = solve_sdp(K, E ,N, M)
%   SOLVE_SDP creates an cvx instance with given dimensions to solve an SDP
%   problem
%   optim_var = solve_sdp(K, E ,N, M) returns a struct with all optimization
%   variables

P_hover = ones([1 N])*1e-3; % Power needed for hovering per time slot
P_max = 5*ones([1 N]); % maximum power of UAV
v = zeros([2 1 N]); % UAV velocity for each time slot
q = zeros([2 1 N]); % UAV position for each time slot
tau = 0.1; % factor for forcing a to be binary
R_min = ones([K 1]); % minium avarage data rate for each user k
Ns_max = 2; % maximum sensing time slots
d_e = [2 -2; 2 2]; % sensing target locations
D = 30; % maximum distance to sensing target
a_t = ones([E N]) * 1/(E*N); % fixed value for taylor approx of a
mu_t = ones([K N]) * 1/(K*N); % fixed value for taylor approx of mu
phi_t = ones([K N]) * 1/(K*N); % fixed value for taylor approx of phi
a_diff = 1;
H = repmat(diag(ones([1 M]) ),[1 1 K N]);

optim_var.W = 0; optim_var.W_t = 0; optim_var.a = 0; optim_var.mu = 0; optim_var.phi = 0;

while (a_diff > 0.1)
    cvx_begin sdp
        variables W(M,M,K,N) W_t(M,M,K,N) a(E,N) mu(K,N) phi(K,N);
        expressions P_t_per_time(1,N) dist_sens(1,N) a_approx sinr_min(K,N) P_rec(K,N) P_int(K,N);
    
        P_t = 0; P_h = 0; P_f = 0;  % variables for transmit, hovering, flying power 
        I_a = 0; % sensing indicators
        
        for n = 1:N
            for k = 1:K % calculation of auxiliary expressions for constraints
                sinr_min(k,n) = 0.5.*( mu(k,n) + phi(k,n) ).^2 - mu_t(k,n) * (mu(k,n) + phi_t(k,n)) ...
                    - phi_t(k,n) * (phi(k,n) - phi_t(k,n));
                P_rec(k,n) = trace(W(:,:,k,n) * H(:,:,k,n));
                interference = 0;
                for m=1:K
                   interference =  interference + trace(W(:,:,m,n) * H(:,:,k,n)); 
                end
                P_int(k,n) = interference - P_rec(k,n);
            end
            P_t_per_time(n) = trace(sum(W(:,:,:,n),3)); % transmitted power each timestep, expression for C1
            dist_sens(n) = sum(a(:,n)' .* vecnorm(q(:,:,n) - d_e)); % distance to sensing target. expression for C7
            I_a = I_a + sum( a(:,n) - 2 .* a_t(:,n) .* ( a(:,n) - a_t(:,n) ) ); % taylor approx, expression for C11b
            
            P_t = P_t + P_t_per_time(n);
            P_h = P_h + sum(a(:,n) * P_hover(n));
            P_f = P_f + (1 - sum( a(:,n) ) ) * calc_p_fly(v(n));
        end
    
        minimize( 1/N * (P_t + P_h + P_f) + tau * I_a)
%           minimize(tau * I_a) 
        	
        subject to
            W_t
%             P_t_per_time <= P_max; % C1
%             P_rec(:) >= sinr_min(:); % C2a
%             P_int(:) <= phi(:); % C2b
%             1/N * sum(log(1 + mu), 2) >= R_min; % C2c
%             sum(a, 1) <= 1; % C5
%             sum(a, 2) <= Ns_max; % C6
%             dist_sens <= D; % C7
%             0 <= a(:) <= 1; % C11a
%             I_a <= 0; % C11b
       
    cvx_end
    
    optim_var.W = W;
    optim_var.W_t = W_t;
    optim_var.a = a;
    optim_var.mu = mu;
    optim_var.phi = phi;

    a_diff = abs(optim_var.a - a_t)
    a_t = optim_var.a;
end
end