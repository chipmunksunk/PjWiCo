K = 2; % communication users
E = 3; % potential sensing targets
N = 10; T = 30; dt = T/N; % Total time T, time slots N, time slot duration dt
M = 5; % number of UAV antennas
w = ones([M K N]); % transmit beamforming vector
W = ones([M M K N]); % Auxiliary variable of w for SDP
a = rand([E N]);
mu = ones([K N]);
d_e = zeros([2 E]);
q = zeros([2 1 N]); % UAV position for each time slot
a_t = ones([E, N]) * 1/(E*N); % fixed value for taylor approx of a
%%
W_e = zeros([M M K E N]);
for e = 1:E
   W_e(:,:,:,e,:) = W;
end
size(W_e)
%%
W_a = W(:,:,:,1);
W_b = sum(W_a,3);
W_c = sum(W,4);
W_d = sum(sum(W,4),3);
%%
a_t = rand([E N]);
a_diff = 1;
cvx_clear
while(a_diff > 1e-3)

    cvx_begin
        variable a(E,N)
        I_a = 0;
        for n=1:N
            I_a = I_a + sum( a(:,n) - a_t(:,n).^2 - 2 .* a_t(:,n) .* ( a(:,n) - a_t(:,n) ) );
        end
        
    %     minimize (sum(sum( a - a.^2, 2), 1))
        minimize(I_a)
    
        subject to
            0 <= a <= 1;  
    cvx_end
    
        a_diff = sum((a-a_t).^2, 'all');
        a_t = a;
end
%%
% close all
% figure('units', 'normalized','OuterPosition',[0 0 1 1])
% set(groot,'defaultfigureunits','normalized','defaultfigureposition',[0 0 1 1])
figure
plot(1:10)


