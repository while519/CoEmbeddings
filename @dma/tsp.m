%From paper: Hahsler et al., Getting Things in Order: An Intro to the R Package seriation, JSS, 2008
%
%D       : an nxn distance matrix.
%'starts': number of generated tours by the heuristic solver tspsearch()
%perm    : the permutation vector

%JY Goulermas, 2013

function [perm, details] = tsp(D, varargin)

  if mod(length(varargin),2)
    error('dma:tsp', 'missing parameters')
  end

  n      = size(D,1);
  starts = n; %default values

  for i = 1 : 2 : length(varargin) - 1
    switch lower( varargin{i} )
      case 'starts'
        starts = min(n, varargin{i+1});
      otherwise
        error('dma:tsp', 'unknown parameter: %s', varargin{i} )
    end
  end

  Solve               = @(x) tspsearch(x, starts+1);             %could use different solvers later
  [perm,details.cost] = Solve([zeros(n+1,1), [zeros(1,n); D] ]); %inserts a dummy node with zero connection costs to all others,
  cut                 = find(perm == 1);                         %to convert the Hamiltonian cycle to a path problem
  perm                = perm([cut+1:end,1:cut-1]) -1;
  
end

%--------------------------------------------------------------------------------------------------

