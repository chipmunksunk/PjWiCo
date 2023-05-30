%% Master file for Project Seminar Wireless Communication
%
%   Implementation of an algorithm for resource allocation and trajectory
%   design for UAV assisted ISAC. Equations are referenced from following paper:
%   https://arxiv.org/abs/2302.10124
%
%% Setup
%
clearvars
close all

K = 2; % communication users
E = 1; % potential sensing targets
N = 2; T = 30; dt = T/N; % Total time T, time slots N, time slot duration dt
M = 2; % number of UAV antennas

optim_var = solve_sdp(K, E ,N, M);

%% check solution
%
W_t = optim_var.W_t;
W = optim_var.W;
a = optim_var.a;
P_max = 1;
W_e = zeros([M M K E N]);

for e = 1:E
    W_e(:,:,:,e,:) = W;
end

for n = 1:N
    for k = 1:K 
        for e = 1:E
            temp1 = - W_t(:,:,k,e,n) + a(e,n) .* P_max .* eye(M); % C12
            try chol(temp1);
                
            catch ME
                disp('Matrix is not symmetric positive definite')
            end
            temp2 = - W_t(:,:,k,e,n) + W_e(:,:,k,e,n); % C13
            try chol(temp2);
                
            catch ME
                disp('Matrix is not symmetric positive definite')
            end
            % C14 is applied when creating variable 
            temp3 = W_t(:,:,k,e,n) - W_e(:,:,k,e,n) + (1-a(e,n)) .* P_max .* eye(M); % C15
            try chol(temp3);
                
            catch ME
                disp('Matrix is not symmetric positive definite')
            end
            temp4 = W_t(:,:,k,e,n);
            try chol(temp4);
                
            catch ME
                disp('Matrix is not symmetric positive definite')
            end

        end
    end
end






