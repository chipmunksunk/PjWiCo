function sv = steer_vec(q, d, M, alti)
% STEER_VEC calculates the steering vector for given UAV position and
% target/user
%
%   q - UAV location
%   d - user location
%   M - # of antennas / length of vector
%   alti - altitude of UAV
if(size(q,1)~=size(d,1))
    ME = MException('MyComponent:noSuchVariable', ...
        'dimensions are mismatched');
    throw(ME)
end

theta = acos(alti/sqrt(norm(q-d)^2 + alti^2));
m = 0:M-1;
sv = exp(1i.*pi.*m.*cos(theta));
end