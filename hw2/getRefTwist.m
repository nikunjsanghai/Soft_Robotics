function refTwist = getRefTwist(d1, tangent, refTwist)
% Compute reference twist
[ne, ~] = size(d1);

% refTwist = zeros(ne+1, 1);

for c=2:ne % all internal nodes
    u0 = d1(c-1,:); % ref director on (c-1)th edge
    u1 = d1(c,:); % ref director on c-th edge
    t0 = tangent(c-1,:);
    t1 = tangent(c,:);

    % Method 1
    % refTwist(c) = computeReferenceTwist(u0, u1, t0, t1);
    
    % Method 2
    ut = parallel_transport(u0, t0, t1);
    ut = rotateAxisAngle(ut, t1, refTwist(c));
    sgnAngle = signedAngle(ut, u1, t1);
    refTwist(c) = refTwist(c) + sgnAngle;
end

end

function vNew = rotateAxisAngle( v, z, theta )
if theta ==0
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v)*(1-c)*z;
end
end
