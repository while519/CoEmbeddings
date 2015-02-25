%QAP-based seriation
%D   : an nxn distance matrix
%perm: the permutation vector

%JY Goulermas, 2013

function [perm, details] = qap(D, varargin)

  n         = size(D,1);
  method    = 'vogelstein'; %default values
  max_iters = 90;
  W         = dma.template(n,'linear');
  
  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'method'
        method = varargin{i+1};
      case 'max_iters'
        max_iters = varargin{i+1};
      case 'w'
        W = varargin{i+1};
      otherwise
        error('dma:qap', 'unknown parameter: %s', varargin{i} )
    end
  end

  switch lower(method)
    case 'vogelstein' %calling qap code used in: Vogelstein et al., Large Graph Matching via Fast Approximate Quadratic Programming, Optimization and Control, 2012
      [details.cost, perm, ~, details.iters,~, ~] = sfw(W, D, max_iters); %minimises sum(sum(W.*D(p,p)))

    case '???'
      perm = 1:n;

    otherwise
      error('dma:qap', 'unknown method: %s', method )
  end
  

end
