%From paper: Hahsler et al., Getting things in order, JSS, 2008
%
%Orders the rows using the first component from the PCA of D (for one-mode),
%and additionally the columns of D (for two-mode) with the first component of D'.
%D   : an nxm matrix with the two-mode data, or a square one-mode matrix.
%mode: 1 or 2 to specify the operating mode
%perm: the permutation vector(s)

%JY Goulermas, 2012

function [perm] = pca(D, mode)

  X        = bsxfun(@minus, D, mean(D,1) );
  [~,~,M]  = svd(X,0);
  [~,perm] = sort( X * M(:,1) );

  if mode == 2
    perm = struct('row', perm, ...
                  'col', dma.pca(D', 1) ...
                 );
  end

end
