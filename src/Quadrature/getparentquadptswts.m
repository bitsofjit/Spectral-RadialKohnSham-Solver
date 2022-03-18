function quadsetup = getparentquadptswts(quadtype, nsp)
% Generate the quadrature points and weights in the biunit parent domain
switch quadtype
    case 1      % Use Gauss quadrature
        xiq = gausspoints(nsp);
        wtq = gaussweights(nsp);
    case 2      % Use Lobatto quadrature
        xiq = lobatto_points(nsp);
        wtq = lobatto_weights(nsp);
end

quadsetup = struct('points',xiq, 'weights',wtq);

return

end        