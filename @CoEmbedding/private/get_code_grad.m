function [f,g] = get_code_grad(x,params)
%
% [f,g] = get_code_grad(x,params)
%
% Calc gradient of the CODE target function.
%

[phi,psis,as,bs] = read_mway_params(x,params);

f = 0;
phigrad =0;
ndist = length(params.pxys);

for i=1:ndist
  [val,gphi,gpsi{i},ga{i},gb{i}] = ...
      call_code_grad_c(phi,psis{i},as{i},bs{i},sparse(params.pxys{i}),...
		       params.px0s{i},params.py0s{i},params.x_cond);
  ga{i} = ga{i}*params.w_nxy(i);
  gb{i} = gb{i}*params.w_nxy(i);

  phigrad = phigrad+params.w_nxy(i)*gphi;
  f = f + params.w_nxy(i)*val;
end

if isfield(params,'pxx') & ~isempty(params.pxx)
  pxx_x0 = sum(params.pxx,2);
  [val,gphi,dummy,tmpga] = ...
      call_code_grad_c(phi,phi,as{1},as{1},...
		       sparse(params.pxx),pxx_x0,pxx_x0,params.x_cond);
   ga{1} = ga{1} + 2*tmpga*params.w_pxx;  
   f = f + params.w_pxx*val;
   phigrad = phigrad + params.w_pxx*gphi;
end


if isfield(params,'pyy') & ~isempty(params.pyy)
   pyy_y0 = sum(params.pyy,1);
  [val,tmp_gpsi,dummy,tmpgb] = ...
      call_code_grad_c(psis{1},psis{1},bs{1},bs{1},...
		       sparse(params.pyy),pyy_y0,pyy_y0,params.x_cond);
   gb{1} = gb{1} + 2*tmpgb*params.w_pyy;  
   f = f + params.w_pyy*val;
   gpsi{1} = gpsi{1} + params.w_pyy*tmp_gpsi;
end

% Construct overall gradient
g  = [phigrad(:)];
if params.b_fix_phi
  g = [];
else
  g = [phigrad(:)];
end

for i=1:length(params.pxys)
  if ~params.b_fix_psi
    g = [g;gpsi{i}(:)];
  end
  
  if params.x_marg=='F'
    g = [g;ga{i}(:)]; 
  end
  
  if  params.y_marg=='F'
    g = [g;gb{i}(:)]; 
  end
end

return
