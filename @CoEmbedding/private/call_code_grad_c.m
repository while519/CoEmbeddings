function [f,gphi,gpsi,ga,gb,logpxy] = call_code_grad_c(phi,psi,a,b,pxy0,px0,py0,x_cond)

phi_tr = phi';
psi_tr = psi';

[i,j,v] = find(pxy0);
psi_exp_d = pxy0*psi;
phi_exp_d = phi_tr*pxy0;

[gphi,gpsi,ga,gb,f] = code_grad(phi,psi,phi_tr,psi_tr,full(px0),full(py0),phi_exp_d,psi_exp_d',a,b,x_cond,i,j,v);

f = -f;


gphi = -gphi';
gpsi = -gpsi';
ga = ga;
gb = gb;

