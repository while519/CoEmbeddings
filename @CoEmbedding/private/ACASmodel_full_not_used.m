
function [X, Y] = ACASmodel_full(R, ratio, p, alpha, beta)

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
  
  R_scale  = sqrt(diag(1./Dx))*R*sqrt(diag(1./Dy));
  [U,S,V]  = svd(R_scale);   
  
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
    if p==1 % when p==1, the largest  singular value is 1 and the corresponding singular vector is constant
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
    if p==1
      U = U(:, 2:1+dim);
      V = V(:, 2:1+dim);
      S = S(2:1+dim, 2:1+dim);
    else
      U = U(:, 1:dim);
      V = V(:, 1:dim);
      S = S(1:dim, 1:dim);
    end
   
  end  
  
  %SS =zeros(size(S));
  %I1 = R_scale*R_scale';
  %I2 = R_scale'*R_scale;
  %for ii=1:length(beta)
   % if ii>2
    %  I1 = I1*I1';
     % I2 = I2*I2';
      %c = max(max(I1(:)), max(I2(:)));
      %SS = SS + beta(ii)*S.^(2^(ii-1))/c;
    %elseif ii==2
     % c = max(max(I1(:)), max(I2(:)));
      %SS = SS + beta(ii)*S.^(2^(ii-1))/c;
   % elseif ii==1
     % c = max(R_scale(:));
     % SS = SS + beta(ii)*S/c;
    %end
 % end
    
  
  X =  diag(1./Dx.^alpha)*U*S.^beta;
  Y =  diag(1./Dy.^alpha)*V*S.^beta;    

end