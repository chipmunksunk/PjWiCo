function sv = steer_vec(q, d, M, alti)
% STEER_VEC calculates the steering vector for given UAV position and
% target/user
%
%   q - UAV location
%   d - user location
%   M - # of antennas / length of vector
%   alti - altitude of UAV

theta = acos(alti/sqrt(norm(q-d)^2 + alti^2));
m = 0:M-1;
sv = exp(1i.*pi.*m.*cos(theta));
end