function [y, X, Y, para, Rqz, Rq] = ACASfun(x, I)
% information onparameters in vector x

% parameter dim 
% parameter p1 for scaling R
% parameter p2 for scaling R
% parameter alpha for scaling eigenvectors 
% parameter beta for scaling vecotrs

  o           = I{1}; 
  paraname    = I{2};
  fixed_para  = I{3};
  marker =1;
  
  for ii = 1 : 2 : length(fixed_para)-1
    switch lower(fixed_para{ii})
      case 'p'
        p = fixed_para{ii+1};
        marker=0;
      case 'p1'
        p1 = fixed_para{ii+1};
      case 'p2'
        p2 = fixed_para{ii+1};
      case 'alpha'
        alpha  = fixed_para{ii+1};
      case 'beta'
        beta = fixed_para{ii+1};
      case 'dim'
        dim = fixed_para{ii+1};
    end
  end


  for ii = 1 : length(paraname)
    switch lower(paraname{ii})
      case 'dim'
        dim = x(ii);
      case 'p'
        p = x(ii);
        marker = 0;
      case 'p1'
        p1 = x(ii);
      case 'p2'
        p2 = x(ii);
      case 'alpha'
        alpha = x(ii);
      case 'beta'
        beta = x(ii);
    end
  end

  if marker==0    
    [X, Y]  = ACASmodel(o.R, o.opt, dim, p, alpha, beta);
  else    
    [X, Y]  = ACASmodel(o.R, o.opt, dim, [p1, p2], alpha, beta);
  end
  
  Rz      = CalRelationXY( X, Y, o.opt_Rz, 'mean');  
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
   
  switch o.obj_option
    case {'quant_integer', 'quant_int', 'quant_adapt',...
          'fast_quant_integer', 'fast_quant_int', 'fast_quant_adapt'}
        Rq      =  ACAS_R(o.R, o.obj_option, o.obj_para);
        [y, Rqz]  = ACAS_obj(Rq, Rz, o.obj_option, o.obj_para);
    otherwise
        [y, Rqz, Rq]  = ACAS_obj(o.R, Rz, o.obj_option, o.obj_para); 
  end
  para.dim =dim;
  if marker ==0
    para.p =p;
  else
    para.p = [p1, p2];
  end
  para.alpha = alpha;
  para.beta = beta;
  para.x = x;
  
end