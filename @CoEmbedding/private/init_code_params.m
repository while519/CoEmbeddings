function [x0,cg_params] = init_code_params(cg_params,phi_init,psi_init,a_init,b_init,noise)
  
embed_dim = cg_params.dim;
x_marg = cg_params.x_marg;
y_marg = cg_params.y_marg;


if isempty(phi_init)
  phi0 = 2*(0.5-rand(cg_params.NX,embed_dim))*noise;
else
  phi0 = phi_init;
end

if ~cg_params.b_fix_phi
  x0 = [phi0(:)];
else
  cg_params.phi0 = phi0;
end

% Prepare marginals
for i=1:length(cg_params.pxys)
  cg_params.pxys{i} = cg_params.pxys{i}/sum(cg_params.pxys{i}(:));
  cg_params.px0s{i} = sum(cg_params.pxys{i},2);
  cg_params.py0s{i} = sum(cg_params.pxys{i},1);      
  cg_params.NYs(i) = length(cg_params.py0s{i});
  cg_params.w_nxy(i) = 1;
  
  if isempty(psi_init)
    psi0 = 2*(0.5-rand(cg_params.NYs(i),embed_dim))*noise; 
  else
    psi0 = psi_init;
  end
  
  if ~cg_params.b_fix_psi
    x0 = [x0;psi0(:)];            
  else
    cg_params.psi0{i} = psi0;
  end
  
  if isempty(a_init)
    a0 = 2*(0.5-rand(cg_params.NX,1))*noise;
  else
    a0 = a_init;
  end
  
  if x_marg=='F'
    x0 = [x0;a0(:)];
  else
    switch x_marg
      case 'U'
	a0 = zeros(cg_params.NX,1);
      case 'M'
	a0 = log(cg_params.px0s{i});
    end
    cg_params.a0{i} = a0;
  end
    
  if isempty(b_init)
    b0 = 2*(0.5-rand(cg_params.NYs(i),1))*noise;
  else
    b0 = b_init{i}
  end
  
  if y_marg=='F'
    x0 = [x0;b0(:)];                
  else
    switch y_marg
      case 'U'
	b0 = zeros(cg_params.NYs(i),1);
      case 'F'
	b0 = log(cg_params.py0s{i});
    end
    cg_params.b0{i} = b0;                
  end
end
