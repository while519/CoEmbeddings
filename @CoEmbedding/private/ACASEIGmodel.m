
function [X, Y] = ACASEIGmodel(R, opt, dim, p, w_between, w_within, alpha)

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
      R_scale  = sqrt(diag(1./Dx))*R*sqrt(diag(1./Dy));
      if w_between+w_within==0
        w_between=1;
        w_within=1;
      end
      T = [w_within*R_scale*R_scale'/ max(max(R_scale*R_scale')), w_between*R_scale/max(R_scale(:)); w_between*R_scale'/max(R_scale(:)), w_within*R_scale'*R_scale/max(max(R_scale'*R_scale))];
    case 'ga'
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
      if w_between+w_within(1) +w_within(2)==0
        w_between=1;
        w_within(1)=1;
        w_within(2)=1;
      end
      T = [w_within(1)*R_scale*R_scale'/ max(max(R_scale*R_scale')), w_between*R_scale/max(R_scale(:)); w_between*R_scale'/max(R_scale(:)), w_within(2)*R_scale'*R_scale/max(max(R_scale'*R_scale))];
    end
  
  
  
  [U,S]  = eig(T); 
  [S0,ind] = sort(real(diag(S)), 1, 'descend');
  U = U(:, ind);
  S = S(ind, ind);
  
  if (length(p)==1 & p(1)==1)
    U = U(:, 2:dim+1);
    S = S(2:dim+1, 2:dim+1);
  elseif p(1)==1&p(2)==1
    U = U(:, 2:dim+1);
    S = S(2:dim+1, 2:dim+1);
  else
    U = U(:, 1:dim);
    S = S(1:dim, 1:dim);
  end
  
  [n,m]  = size(R);
  V = U*sqrt(S);  
  X = real(diag(1./Dx.^alpha)*V(1:n,:));
  Y = real(diag(1./Dy.^alpha)*V(n+1:end,:));
  

end