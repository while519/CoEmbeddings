%Creates a flow matrix W for QAP minimisation with a dissimilarity matrix D, with the purpose of seriating D.
%The sign of W is adjusted to keep it a minimisation problem
%n   : side length of the matrix
%type: string of the name for the type of the matrix
%k   : additional parameter: for 'band' k>0 is the band radius

%JY Goulermas & A Kostopoulos, 2012

function [W] = template(n, type, k)

  if ~exist('type', 'var')
    type = 'linear';
  end
  
  [C,R]= meshgrid(1:n);
  
  switch lower(type)
  
    case 'linear'%-least squares
      W = -abs(R-C);

    case 'squared' %squared (inertia)
      W = -(R-C).^2;

    case 'circular'
      t      = (n - mod(n,2)) / 2;
      W      = abs(R-C);
      W(W>t) = n-mod(n,2) - W(W>t);
      W      = -W;
      
    case 'band'
      W      = abs(R-C);
      W(W>k) = 0;
      W      = -W;

    case 'tsp'
      W = abs(R-C)==1;
      
    case 'rank1'
      u = linspace(-1, +1, n)';
      W = u*u';

    case 6%???
      W = double('?');

    otherwise
      error('unknown seriation matrix type')
  end
  
 W = double(W);

end