function refTwist = computeReferenceTwist(u1, u2, t1, t2)
% refTwist is the reference twist to move u1 with tangent t1 to u2 with
% tangent t2.
ut = parallel_transport(u1, t1, t2);
refTwist = signedAngle(ut, u2, t2);
end
