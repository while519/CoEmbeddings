
function [X, Y, svdres] = ACASmodel(R, opt, ratio, p, alpha, beta, svdres)

  if nargin==6
    svdres=[];
  end

  if isempty(svdres)
    switch opt
      case {'grid_search', 'grid'}  
        if p == 0
          Dx=1;
          Dy=1;
        elseif p == -1
          Dx = max(R');
          Dy = max(R );
        else
          Dx = sum(R.^p,2).^(1/p);
          Dy = sum(R.^p,1).^(1/p);
        end
      case 'ga'
        switch length(p)
          case 2
            if p(1) == 0
              Dx=1;
            elseif p(1)<0
              Dx = max(R');
            else
              Dx = sum(R.^p(1),2).^(1/p(1));
            end

            if p(2) == 0
              Dy=1;
            elseif p(2)<0
              Dy = max(R );
            else
              Dy = sum(R.^p(2),1).^(1/p(2));
            end
          case 1
            if p == 0
              Dx=1;
              Dy=1;
            elseif p<0
              Dx = max(R');
              Dy = max(R );
            else
              Dx = sum(R.^p,2).^(1/p);
              Dy = sum(R.^p,1).^(1/p);
            end            
        end
    end
    R_scale  = sqrt(diag(1./Dx))*R*sqrt(diag(1./Dy));    
    [U,S,V]  = svd(R_scale);   
    svdres = {U,S,V, Dx, Dy};
  else
    U = svdres{1};
    S = svdres{2};
    V = svdres{3};
    Dx = svdres{4};
    Dy = svdres{5};
  end
  
  
  
  [n,d] = size(R);  
  L = diag(S) .^2 / (n-1);
  total_var = sum(L);
  if strcmpi(ratio,'full') %put even the zero elements when rank-deficient covariance
    rank_max = min(n-1,d);
    if rank_max  < d
       U(:, rank_max+1 : end) = 0;
       V(:, rank_max+1 : end) = 0;
       S(rank_max+1 : end, rank_max+1 : end) = 0;
    end
    
    if (length(p)==1 & p(1)==1)
      U = U(:, 2:end);
      V = V(:, 2:end);
      S = S(2:end, 2:end);
    end
    
    if p(1)==1&p(2)==1
      U = U(:, 2:end);
      V = V(:, 2:end);
      S = S(2:end, 2:end);
    end
    
  else
    if strcmpi(ratio,'max') %all the effective elements only
      dim = rank( cov(R) );
    elseif ratio > 0.0 & ratio < 1.0 %only elements capturing the significant variance
      var    = L / sum(L);
      hist   = cumsum(var); %cumulative histogram of increasing eigenvalues
      dim = sum(hist <= ratio); %find indices of significant components and how many
      dim = max(dim, 1); %at least one axis
    elseif ratio >= 1 & ceil(ratio) == ratio %exact number of components needed
      dim = ratio;
    else
      error('incorrect variance ratio');
    end
        
    if (length(p)==1 & p(1)==1)
      tt = min(1+dim, size(U,2));
      tt = min(tt, size(V,2));
      U = U(:, 2:tt);
      V = V(:, 2:tt);
      S = S(2:tt, 2:tt);
    elseif p(1)==1&p(2)==1
      tt = min(1+dim, size(U,2));
      tt = min(tt, size(V,2));
      U = U(:, 2:tt);
      V = V(:, 2:tt);
      S = S(2:tt, 2:tt);
    else
      tt = min(dim, size(U,2));
      tt = min(tt, size(V,2));
      U = U(:, 1:tt);
      V = V(:, 1:tt);
      S = S(1:tt, 1:tt);
    end    
  end  
  
  X =  real(diag(1./Dx.^alpha)*U*S.^beta);
  Y =  real(diag(1./Dy.^alpha)*V*S.^beta);    


end