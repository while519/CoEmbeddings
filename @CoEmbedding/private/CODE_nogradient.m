function [X, Y] = CODE(R, dim, ini)
  pxys =R;
  embed_dim= dim;
  
  if nargin<3

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

    [x0,cg_params] = init_code_params(cg_params,phi_x_init, ...
                psi_y_init,a_init,b_init,noise);
    x0= x0(1:(size(R,1)+size(R,2))*dim);         
  else
    x0= [ini{1}(:); ini{2}(:)];
  end
  
  o.R = R;
  o.dim =dim;   
  o.rowa =  size(R,1);
  o.rowb =  size(R,2);
  o.mark = size(R,1)*dim;
  
  for ii=1:length(x0);
    range{ii} = [0,1];
  end
  
  setting = {'Function2Optimize',     @myfunc,...
             'FixedArgument',         o,...
             'ContinuousVariable',    range,...
             'PopulationSize',        30,...
             'MaxGen',                500,...
             'IniPopulation',         x0',...
             'toel',                  1e-10,...
             'display',               1,...
             };
           
  obj = GAParaOpt(setting);
  obj = optimize(obj);
  res  = obj.output;
  x=res. Optimal_Cpara; 
  
  %x = [];
  %for ii=1:10
   % if ~isempty(x)
    %   x0 = x;
    %end
    %[x,fval,exitflag] = fminunc(@(x)myfunc(x, o),x0,options);
  %end

  a = x(1:o.mark); b=x(o.mark+1:end);
  X = vec2mat_tt(a, o.rowa);
  Y = vec2mat_tt(b, o.rowb);

end



%==========
function y = myfunc(x, o)
  a = x(1:o.mark); b=x(o.mark+1:end);
  X = vec2mat_tt(a, o.rowa);
  Y = vec2mat_tt(b, o.rowb);
  y = Log_Likelihood(X, Y, o.R);
end

%=============
function X = vec2mat_tt(x, row)
  n   = length(x);
  col = n/row;
  for ii=1:col
    X(1:row,ii) = x(1+(ii-1)*row:row+(ii-1)*row);
  end
end
%==============
function o = Log_Likelihood(X, Y, R)

  R = R/sum(R(:));

  Px = sum(R,2);
  Py = sum(R,1);
  P = CalRelationXY(X, Y, 'gaussian2',1);
  D = CalRelationXY(X, Y, 'sqeuclidean');
  Zx = sum(P*diag(Py),2);
  
  item1 = R.*D; item1 = sum(item1(:));
  item2 = sum(log(Zx).*Px);
 
  o = item1 +item2;
  
  %Dx = diag(sum(R,2));
  %Dy = diag(sum(R,1));
  %P = CalRelationXY(X, Y, 'gaussian2',1);
  %Zx = diag(1./sum(P*Dy,2));
  %Pcm = Zx*Dx*P*Dy;

  %o = -R.*log(Pcm);
  %o = sum(o(:));

end


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
end
% ===============================================================
% Check some arguments
function check_args(phi_x_init)
if (iscell(phi_x_init)), error('phi_x_init cannot be a cell array');end
end
