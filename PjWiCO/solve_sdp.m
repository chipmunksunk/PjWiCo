function optim_var = solve_sdp(K, E ,N, M)
%   SOLVE_SDP creates an cvx instance with given dimensions to solve an SDP
%   problem
%   optim_var = solve_sdp(K, E ,N, M) returns a struct with all optimization
%   variables
%
%   K - # of comm users
%   E - # of sensing targets
%   N - # of time slots
%   M - # of antennas

P_hover = ones([1 N])*1e-3; % Power needed for hovering per time slot
P_max = 10; % maximum power of UAV in Watts
v = zeros([N 2]); % UAV velocity for each time slot
q = zeros([N 2]); % UAV position for each time slot
altitude = 40; 
sigma_e2 = 1e-8; % target noise power in Watts
sigma_k2 = sigma_e2;
beta_o = 1e-3;
tau = P_max*10; % factor for forcing a to be binary
R_min = 1e-2; % minium avarage data rate
snr_echo_min = 0; % minimum comm sinr of 0 dB
Ns_max = 1; % maximum sensing time slots
d_e = [-5 5; 10 10]; % sensing target locations
d_k = [0 5; 20 20; -10 10]; % com target locations
D = 50; % maximum distance to sensing target
a_t = ones([E N]); % fixed value for taylor approx of a
mu_t = zeros([K N]) * 1/(K*N); % fixed value for taylor approx of mu
phi_t = zeros([K N]) * 1/(K*N); % fixed value for taylor approx of phi
a_diff = 1; mu_diff = 1; phi_diff = 1; % variables for SCA of taylor approx
H = repmat(eye(M),[1 1 K N]); % channel matrix
for n =1:N
    for k = 1:K
        hk = steer_vec(q(n,:), d_k(k,:), M, altitude) * beta_o / sqrt( norm( q(n,:) - d_k(k,:) )^2 + altitude^2 );
        H(:,:,k,n) = hk'*hk;
        H(:,:,k,n) = H(:,:,k,n) * 1e+9;
    end
end


optim_var.W = 0; optim_var.W_t = 0; optim_var.a = 0; optim_var.mu = 0; optim_var.phi = 0;
cvx_clear

while (a_diff > 1e-3 || mu_diff > 1e-3 || phi_diff > 1e-3)
    cvx_begin sdp
        variable W(M,M,K,N) semidefinite;
        variable W_t(M,M,K,E,N) semidefinite;
        variables a(E,N) mu(K,N) phi(K,N);
        expressions P_t_per_time(1,N) dist_sens(1,N) nu(K,N) P_rec(K,N) P_int(K,N) snr_echo(E,1);

        P_t = 0; P_h = 0; P_f = 0;  % variables for transmit, hovering, flying power 
        I_a = 0; % sensing indicators
        
        for n = 1:N
            for k = 1:K % calculation of auxiliary expressions for constraints

                P_rec(k,n) = trace(W(:,:,k,n) * H(:,:,k,n)); % received power for each user k and time slot n, expr. for C2a
                nu(k,n) = 0.5 * ( mu(k,n) + phi(k,n) )^2 - 0.5 * ( mu_t(k,n)^2 + phi_t(k,n)^2 )...
                    - mu_t(k,n) * ( mu(k,n) - mu_t(k,n) )...
                    - phi_t(k,n) * ( phi(k,n) - phi_t(k,n) ); % minimum P_rec for each user k and time slot n, expr. for C2a
                int_temp = 0;
                for j = 1:K
                    int_temp = int_temp + trace(W(:,:,j,n) * H(:,:,k,n));
                end
                P_int(k,n) = int_temp - P_rec(k,n); % interference power at user k, expr. for C2b
            end
            P_t_per_time(n) = trace(sum(W(:,:,:,n),3)); % transmitted power each timestep, expression for C1
            I_a = I_a + sum( a(:,n) - a_t(:,n).^2 - 2 .* a_t(:,n) .* ( a(:,n) - a_t(:,n) ) ); % taylor approx, expression for C11b
            
            P_t = P_t + P_t_per_time(n);
            P_h = P_h + sum(a(:,n) * P_hover(n));
            P_f = P_f + (1 - sum( a(:,n) ) ) * calc_p_fly(v(:,n));
        end
        
        
        for e = 1:E
          snr_echo_temp = 0;
            for n = 1:N
                sv = steer_vec(q(n,:), d_e(e,:), M, altitude);
                % snr_echo_temp = snr_echo_temp + beta_o^2 * sv * sum(W_t(:,:,:,e,n), 3) * sv' / ( 16*pi*sigma_e2*sqrt( norm(q(n,:)-d_e(e,:))^2 + altitude^2 )^4 );
                snr_echo_temp = snr_echo_temp + sv * sum(W_t(:,:,:,e,n), 3) * sv' / ( sigma_e2 * sqrt( norm(q(n,:)-d_e(e,:))^2 + altitude^2 )^4 );
            end
            snr_echo(e) = snr_echo_temp;
        end
        
    
        minimize( 1/N * (P_t + P_h + P_f) + tau * I_a)
%         minimize(I_a)
        %%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%	
        subject to

        for n = 1:N
            for k = 1:K 
                for e = 1:E
                    W_t(:,:,k,e,n) <= a(e,n) .* P_max .* eye(M); % C12
                    W_t(:,:,k,e,n) <= W(:,:,k,n); % C13
                    % C14 is applied when creating variable 
                    W_t(:,:,k,e,n) >= W(:,:,k,n) - (1-a(e,n)) .* P_max .* eye(M); % C15
                end
            end
        end 
            % P_t_per_time <= P_max; % C1 
            P_rec(:) >= nu(:); % C2a
            P_int(:) <= phi(:); % C2b
            1/N * sum_log(1 + mu, 2) / log(2) >= R_min*ones([K 1]); % C2c
            snr_echo >= snr_echo_min; % C4
            % sum(a, 1) <= 1; % C5
            % sum(a, 2) <= Ns_max; % C6
            % 0 <= a(:) <= 1; % C11a
        
            
       
    cvx_end
    
    optim_var.W = W; optim_var.W_t = W_t; optim_var.a = a; optim_var.mu = mu; optim_var.phi = phi;

    a_diff = sum((optim_var.a - a_t).^2,'all')
    mu_diff = sum((optim_var.mu - mu_t).^2, 'all')
    phi_diff = sum((optim_var.phi - phi_t).^2, 'all')
    a_t = optim_var.a; mu_t = optim_var.mu; phi_t = optim_var.phi;

end
end