function [phi_x,psi_y,out] = code(pxys,embed_dim,options)
%
% function [PHI_X,PSI_Y]=CODE(PXY,EMBED_DIM,OPTIONS)
%
% Implement Euclidean embedding of co-occurence data, using
% conjugate gradient.
%
% Inputs: 
% =======
% 
% PXY       - Joint probability matrix. Can also be a cell array of
%             co-occurrence joint probability matrices. In this case
%             all matrices should have a common X (same number of rows).
%
% EMBED_DIM - Dimension of the low-dimensional embedding.
% 
% OPTIONS   - a structure with the following fields.
%      
%     All options have defaults. 
% 
%     n_restarts   - repeat the algorithm n_restarts times and choose
%                    the one that maximizes the likelihood
%                    (default=20);
%     conjgrad_opts- options vector for conjgrad. See CONJGRAD
%                    (defaults described in CONJGRAD help).
%     w_nxy        - weight of each co-occurrence matrix. Only
%                    relevant for the multiple inputs
%                    case. Determines the relative weight of each
%                    joint distribution pair in the likelihood
%                    maximization. (default=Uniform)
%     prt_level    - verbosity. (default=1);
%     seeds        - vector of seeds for the random generator, a
%                    seed for each restart. default = [1:n_restarts];
%     
%
% Outputs: 
% ========
% phi_x - Embeding of x objects. A matrix of size (|X|,embed_dim)
% psi_y - Embeding of y objects. If pxy is a cell array, psi_y will
%         be a cell array of the same size. Otherwise, it will be
%         the matrix of y embedings of size (|Y|,embed_dim).
% out (optional)
%         A structure with the following fields: 
%    as  - the optimal As vector
%    bs  - the optimal As vector
%    lik - the Likelihood by the model
% 
% See
% Globerson A, Chechik G, Pereira F and Tishby N,  
% Euclidean embedding of co-occurence data, 
% Adavances in Neural Information Processing systems 17, (NIPS*2004)
% 
% (C) G. Chechik, A. Globerson 2004
% Comercial usage requires written permission.
% Used conjgrad.m and linesearch.m by Hans Bruun Nielsen with minor
% modifications

% Additional options. These were not covered in our NIPS paper, but
% provide interface to extended models, including non-conditional
% models with several methods for handling the marginals.
% The model considered here are of the form 
%      p(x,y)= 1/Z exp(-|\phi(x)-\psi(y)|^2+A(x)+B(y))
% They differ in the initialization of A(x),B(y) and whether these
% variables are optimized over or fixed.
% 
%     x_cond       - Make model conditional on x (default=1)
%     x_marg       - How to treat marginal of x. only relevant in
%                    the case of non-conditional models (c_cond=0)
%                    One of:  
%                    'F' - free. Optimized over
%                    'U' - fixed as Uniform
%                    'M' - fixed as the empirical marginals
%     y_marg       - How to treat marginal of y. see x_marg.
%     a_init       - Value for initializing a(x) to. Will only be 
%                    used if x_marg='F'
%     b_init       - Value for initializing b(y) to. Will only be 
%                    used if y_marg='F'
%     pxx          - It is possible to add a distribution p(x,x)
%                    which specifies co-occurence of x's. It will
%                    be modeled as 
%                    p(x1,x2)=1/Z exp(-|\phi(x1)-\phi(x2)|^2+A(x1)+A(x2))
%     w_pxx        - The weight assigned to the likelihood of the
%                    p(x,x) model in the likelihood maximization
%     pyy          - See pxx
%     w_pyy        - See w_pxx
%     phi_x_init   - A fixed value to set phi(x) to. If phi_x_init does
%                    not exist or is empty, phi(x) is optimized
%                    over (the default).
%     b_fix_phi    - Says if we should optimize over phi or keep it
%                    fixed. 
%     psi_y_init   - Value to initialize psi(y) to. 
%     b_fix_psi    - Says if we should optimize over psi or keep it
%                    fixed. 
%     noise        - Parameters will be randomly initialized to a
%                    uniform variable between -noise and +noise (default=1)


% Handle the case of a single joint distribution
if ~iscell(pxys)
   npxys = {};
   npxys{1} = pxys;
   clear pxys;
   pxys = npxys;
   clear npxys;
end

if nargin<3
  options=[];
end

% Parse user options
[options,n_restarts]       = take_from_struct(options,'n_restarts',10);
[options,prt_level]        = take_from_struct(options,'prt_level',1);
[options,cg_params.w_nxy]  = take_from_struct(options,'w_nxy',[]);
[options,cg_opts]          = take_from_struct(options,'conjgrad_opts',-1*ones(1,9));
[options,seeds]         = take_from_struct(options,'seeds',[1:n_restarts]);
seeds(end+1:n_restarts) = [length(seeds)+1:n_restarts];

% Parse undocumented options 
[options,phi_x_init]    = take_from_struct(options,'phi_x_init','');
[options,x_marg]        = take_from_struct(options,'x_marg','F');
[options,y_marg]        = take_from_struct(options,'y_marg','F');
[options,noise]         = take_from_struct(options,'noise',1);
[options,psi_y_init]    = take_from_struct(options,'psi_y_init','');
[options,a_init]        = take_from_struct(options,'a_init','');
[options,b_init]        = take_from_struct(options,'b_init','');
[options,pxx]           = take_from_struct(options,'pxx','');
if ~isempty(pxx);  pxx  = pxx/sum(pxx(:)); end
[options,pyy]           = take_from_struct(options,'pyy','');
if ~isempty(pyy);  pyy  = pyy/sum(pyy(:)); end
[options,w_pxx]         = take_from_struct(options,'w_pxx',0);
[options,w_pyy]         = take_from_struct(options,'w_pyy',0);
[options,save_tmp]      = take_from_struct(options,'save_tmp',0);
[options,cg_params.b_fix_phi] = take_from_struct(options,'b_fix_phi',0);
[options,cg_params.b_fix_psi] = take_from_struct(options,'b_fix_psi',0);
[options,cg_params.x_cond]    = take_from_struct(options,'x_cond',0);

% Check some arguments
check_args(phi_x_init);

% Handle conditional model
if (options.x_cond)
  x_marg = 'U';
  y_marg = 'M';
end

% Set cg_params 
[NX,NY] = size(pxys{1});
cg_params.NX = NX;
cg_params.NYs = NY;
cg_params.dim = embed_dim;
cg_params.pxys = pxys;
cg_params.pxx = pxx;
cg_params.pyy = pyy;
cg_params.w_pxx = w_pxx;
cg_params.w_pyy = w_pyy;
cg_params.x_marg = x_marg;
cg_params.y_marg = y_marg;

% Loop over repeats 
for(i_repeat=1:n_restarts)
  seed = seeds(i_repeat);
  % rand('state',seed);
  if(prt_level)
    fprintf('i_repeat = %3d(%d) seed=%d\n', i_repeat,n_restarts,seed);
  end
  [x0,cg_params] = init_code_params(cg_params,phi_x_init, ...
				    psi_y_init,a_init,b_init,noise);
  
  x = conj_grad('get_code_grad',cg_params,x0,cg_opts);
  [curr_phis,curr_psis,curr_as,curr_bs] = read_mway_params(x,cg_params);
  
  likelyhood =  -get_code_grad(x,cg_params);
  if(prt_level)fprintf('\tLikelihood=%10.8f\n',likelyhood); end
    
  % Save best repeat
  liks(i_repeat) = likelyhood;    
  if likelyhood>=max(liks)
    fprintf('\tFound max %g\n',liks(i_repeat));
    bst_phi = curr_phis;
    bst_psi = curr_psis;
    bst_as = curr_as;
    bst_bs = curr_bs;
  end
    
  % Save-to-gile current repeat
  if(save_tmp)
    save_filename = sprintf('tmp_code_seed%d',seed);
    save(save_filename,'curr_phis', 'curr_psis','curr_as', ...
	 'curr_bs','likelyhood')
    fprintf('\tResults for seed %d were saved to "%s.mat"\n', ...
	    seed,save_filename);
  end
  
  
end

liks = full(liks);
[max_lik,bst_repeat] = max(liks);
bst_repeat = bst_repeat(1);
fprintf('Choose best score out of repeats\n');  
fprintf('max_lik=%f bst_repeat=%d\n',max_lik,bst_repeat);    


lik = max_lik;
phi_x = bst_phi;

% Handle single input case 
if length(pxys)==1
  psi_y = bst_psi{1};
else
  psi_y = bst_psi;
end

as = bst_as;
bs = bst_bs;

% Set outputs
out.as = as;
out.bs = bs;
out.lik = lik;

return
  
% ===============================================================
function [out_options,val]=take_from_struct(options,fieldname,default)
%
% Take values from the options structure
%
out_options = options;    
if(isfield(options,fieldname));
    val = getfield(options,fieldname);
else
    val=default;
    out_options=setfield(out_options,fieldname,val);  
end

return

% ===============================================================
% Check some arguments
function check_args(phi_x_init)
if (iscell(phi_x_init)), error('phi_x_init cannot be a cell array');end
return
