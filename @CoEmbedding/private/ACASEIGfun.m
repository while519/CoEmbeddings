function [y, X, Y, para, Rqz, Rq] = ACASEIGfun(x, I)
% information onparameters in vector x

% parameter dim 
% parameter p1 for scaling R
% parameter p2 for scaling R
% parameter alpha for scaling eigenvectors 
%
% parameter w_between for scaling  between-object similarities
%             0             w_between*R
%           w_between*R'         0
% parameter w_within for scaling  within-object similarities
%         w_within(1)*RR'      0
%             0      w_within(2)*R'R


  o           = I{1}; 
  paraname    = I{2};
  fixed_para  = I{3};
  
  
  for ii = 1 : 2 : length(fixed_para)-1
    switch lower(fixed_para{ii})
      case 'p1'
        p1 = fixed_para{ii+1};
      case 'p2'
        p2 = fixed_para{ii+1};
      case 'alpha'
        alpha  = fixed_para{ii+1};
      case 'w_between'
        w_between = fixed_para{ii+1};
      case 'w_within1'
        w_within1 = fixed_para{ii+1};
      case 'w_within2'
        w_within2 = fixed_para{ii+1};
      case 'dim'
        dim = fixed_para{ii+1};
    end
  end

  
  for ii = 1 : length(paraname)
    switch lower(paraname{ii})
      case 'dim'
        dim = x(ii);
      case 'p1'
        p1 = x(ii);
      case 'p2'
        p2 = x(ii);
      case 'alpha'
        alpha = x(ii);
      case 'w_between'
        w_between = x(ii);
      case 'w_within1'
        w_within1 = x(ii);
      case 'w_within2'
        w_within2 = x(ii);
    end
  end
  
  
  [X, Y]  = ACASEIGmodel(o.R, o.opt, dim, [p1, p2], w_between, [w_within1, w_within2], alpha);
  Rz      = CalRelationXY( X, Y, 'gaussian', 'mean');  
  if sum( o.R(:))==1
    Rz = Rz/(sum(Rz(:)));
  end
  if mean(sum(o.R,2))==1
    Dx = sum(Rz,2);
    Rz = diag(1./Dx)*Rz;
  end
  if mean(sum(o.R))==1
    Dy = sum(Rz);
    Rz = Rz*diag(1./Dy);
  end              
   
  [y, Rqz, Rq ]   =  ACAS_obj(o.R, Rz, o.obj_option, o.obj_para);   
  
  para.dim =dim;
  para.p = [p1, p2];
  para.alpha = alpha;
  para.w_between = w_between;
  para.w_within =[w_within1, w_within2];
  para.x = x;
  
end