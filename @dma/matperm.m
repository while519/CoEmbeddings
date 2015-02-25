%Converts a permutation from a vector/matrix to a matrix/vector representation.
%No input validation is performed.

%JY Goulermas & A Kostopoulos, 2013

function [out] = matperm(in)

  [n,b] = size(in);

  if isvector(in)
    out = sparse(1:max(n,b), in, 1);
  else
    [out,~] = find(in');
  end

end

