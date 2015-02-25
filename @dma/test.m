%From paper:
%
%D   : a nxn dissimilarity matrix
%type: 'sts' or 'neighborhood'
%perm: the permutation vector

%JY Goulermas, 2011

function [perm, details] = spin(D, type)

  if nargin == 1
    type = 'sts';
  end
  
  details = [];

  type = strmatch( lower(type), strvcat('sts', 'neighbourhood', 'fast_neighbourhood', 'rank1', 'rank2') ); if isempty(type), type=0; end
  switch type
    case 1
      [perm, details.iters] = SideToSide(D);
    case {2,3}
      [perm, details.iters] = Neighbourhood(D, type);
    case 4
      [perm] = Rank1(D);
    case 5
      [perm] = Rank2(D);
    otherwise
      error('unknown spin version')
  end
  
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Neighbourhood version
function [perm, iter] = Neighbourhood(D, type)

  n         = size(D,1);
  [r,c]     = meshgrid(1:n);
  sigma     = 25;
  W0        = -(r-c).^2 / n;
  resol     = min(D(:)) /n; %to speedup the LAP (sometimes)
  P         = speye(n);
  cost_prev = inf;
  stop      = false;
  iter      = 0;
  while ~stop
    iter      = iter + 1;
    W         = exp(W0 / sigma);
    [~, ~, W] = SinkhornKnopp(W); %normalise to doubly stochastic (improves mainly the 'fast' case)
    sigma     = max(0.2, sigma / 2); %reduce sigma gradually
    A         = D * P' * W;
    switch type
      case 2 %solve the LAP with p being the col permutations, trace(A*P) == trace(A(:,p)) == cost
        [p,cost] = lapjv(A, resol);
      case 3 %approximate LAP
        [~,q]  = min(A, [], 2); %get minimum of each row and then break ties
        cur    = 1;
        labels = unique(q);
        p      = NaN(n,1);
        for k = labels'
          pos    = find(q==k);
          p(pos) = [cur : cur+numel(pos)-1]';
          cur    = cur + numel(pos);
        end%for
        cost = trace( A(:,p) );
    end%switch
    P         = sparse(p, 1:n, 1); %convert col perm vec to col perm matrix; this still minimises trace(P*A) (not really need to convert between P & p, but kept here for convenience)
    stop      = abs(cost-cost_prev) < 1e-7;
    %imagesc(P*D*P'); abs(cost-cost_prev), pause %%%debug
    cost_prev = cost;
  end
  [perm,~,~] = find(P'); %converts row perm matrix to row perm vec
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Side-to-side version
function [perm, iter] = SideToSide(D)
  n      = size(D,1);
  w      = [1:n]' - (n+1)/2; %implying W=w*w'
  P_prev = speye(n);
  stop   = false;
  iter   = 0;
  while ~stop
    iter   = iter + 1;
    s      = D * P_prev' * w;
    [~, p] = sort(s, 'descend');
    P      = sparse(1:n, p, 1);
    stop   = norm(P*s - P_prev*s, inf) < 1e-7;
    %imagesc(P*D*P'), obj = @(T) w'*T*D*T'*w; obj(P)-obj(P_prev), pause %%%debug
    P_prev = P;
  end
  [perm,~,~] = find(P'); %convert row perm. matrix to perm. vector

  %cost = w'*P*D*P'*w
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Rank-1 approximation (my version)
function [perm] = Rank1(D)
  n     = size(D,1);
  w     = [1:n]' - (n+1)/2;
  [M,L] = EigSorted(D, 'ascend');
  
  switch 002
    case 1 %just get the most optimal (in the relaxed sense)
      z     = M(:,1);
      [~,p] = sort(z);

    case 2 %just get the most similar to w
      for i = 1 : 20
        dev(i) = norm( sort(M(:,i)) - w );
      end
      [~,i] = min(dev);
      z     = M(:,i);
      [~,p] = sort(z);
      
    case 3 %combine the first k??????? - experimental
      k = [3,4];
      MM = M(:,k);
      [MM, I] = sort(MM);
      II = mean(I,2);
      [~,p] = sort(II);
     rp(p)=[1:n];p=rp';
      

%       c    = pinv(MM) * w;
%       z    = I * c;
%      [~,p] = sort(z);
%      
%      II=mean(I,2); [~,p] = sort(II);
%      rp(p)=[1:n];p=rp';
  end
  
  perm = p;
   
  
  P    = sparse(1:n, perm, 1);
  cost = w'*P*D*P'*w
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Rank-2 approximation (our version)
function [perm] = Rank2(D)
  n       = size(D,1);
  [r,c]   = meshgrid(1:n);
  W       = -(r-c).^2 / n;
  %W       = exp(W/88);
  [U,S,V] = svd(W); 
  
  w2      = U(:,1); s2 = U(1,1)./V(1,1) * S(1,1);
  w1      = U(:,2); s1 = U(1,2)./V(1,2) * S(2,2);
  w3      = U(:,3); s3 = U(1,3)./V(1,3) * S(3,3);

  %imagesc(W)
  %imagesc(s1*w1*w1'+s2*w2*w2'), axis square
  %imagesc(s1*w1*w1'+s2*w2*w2'+s3*w3*w3'), axis square

%DO MESS
[w p2] = sort(w2);
  
  
  P_prev = speye(n);
  stop   = false;
  iter   = 0;
  while ~stop
    iter   = iter + 1;
    s      = D * P_prev' * w;
    [~, p] = sort(s, 'descend');
    P      = sparse(1:n, p, 1);
    stop   = norm(P*s - P_prev*s, inf) < 1e-7;
    %imagesc(P*D*P'), obj = @(T) w'*T*D*T'*w; obj(P)-obj(P_prev), pause %%%debug
    P_prev = P;
  end
  

%UNDO MESS
P2 = sparse(1:n, p2, 1);
P = P2'*P;


  [perm,~,~] = find(P'); %convert row perm. matrix to perm. vector

  %cost = w'*P*D*P'*w


end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~












































