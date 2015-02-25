function [X_best, Y_best, time, para] = ACASeig(R, dim, options, matrix_match, ini,display)
% For ga optimisation, the dim, p, alpha of ACAS is either defined as a fixed number or a
% range. For beta, it is either defined as [] (to be searched from [0,1]) or fixed vector with length
% equal to sim_order.
% For grid search optimisation, the dim, p, alpha, beta of ACAS is either defined as a fixed number or a
% set of candidate parameters, because of grid search.

  p1=[];
  p2=[];
  w_within1=[];
  w_within2=[];
  obj_para = [];
  obj_option = 'quant_integer';
  opt = 'ga';
  for ii = 1 : 2 : length(matrix_match)-1
    switch lower(matrix_match{ii})
      case 'obj_para'
        obj_para = matrix_match{ii+1};
      case 'obj_option'
        obj_option = matrix_match{ii+1};
      case {'optimization', 'optimisation'}
        opt = matrix_match{ii+1};
      otherwise
        error( [ 'unknown setting for matrix matching ' matrix_match{ii} ] );
    end
  end
  
  %Default parameter setting
  switch lower(opt)
    case {'grid_search', 'grid'}  
      p               = [-1, 0, 1, 2, 3];
      alpha           = 0:0.5:3;
      w_between       = 0:0.2:1;
      w_within        = 0:0.2:1;  
      if isempty(dim)
        temp = min(size(R,1), size(R,2));
        dim = 2: floor((temp-2)/20)+1 : temp ;
      end
    case 'ga'
      p           = [-1, 10];
      alpha       = [0, 3];
      w_between   = [0, 1];
      w_within    = [0, 1];
      if isempty(dim)
        dim = [2, min(size(R,1), size(R,2))];
      end
  end

  for ii = 1 : 2 : length(options)-1
    switch lower(options{ii})
      case 'p'
        p = options{ii+1};        
      case 'p1'
        p1 = options{ii+1};        
      case 'p2'
        p2 = options{ii+1};
      case 'w_between'
        w_between  = options{ii+1};
      case 'w_within'
        w_within = options{ii+1};        
      case 'w_within1'
        w_within1 = options{ii+1};        
      case 'w_within2'
        w_within2 = options{ii+1};
      case 'alpha'
        alpha = options{ii+1};
      otherwise
        error( [ 'unknown model variables: ' options{ii} ] );
    end
  end
  
  fprintf('ACAS Model Optimisation:\n')
  time = cputime;  
  switch  lower(opt)
    case {'grid_search', 'grid'}     
      best = 1e10;
      perf = [];     
      
      for ii=1:length(p)
        for jj=1:length(w_between)
          for kk=1:length(w_within)  
            for ss=1:length(alpha)
              for tt =1:length(dim)
              
                [X, Y] = ACASEIGmodel(R,  opt, dim(tt),  p(ii),  w_between(jj),  w_within(kk), alpha(ss));              
                Rz     = CalRelationXY(X, Y, 'gaussian', 'mean');                
                if sum( R(:))==1
                  Rz = Rz/(sum(Rz(:)));
                end
                if mean(sum(R,2))==1
                  Dx = sum(Rz,2);
                  Rz = diag(1./Dx)*Rz;
                end
                if mean(sum(R))==1
                  Dy = sum(Rz);
                  Rz = Rz*diag(1./Dy);
                end
                [score, Rq, Rqz]   = ACAS_obj(R, Rz, obj_option, obj_para);  
                if score  < best
                  best       = score ;
                  ind        = [ii,jj,kk,ss,tt];
                  X_best     = X;
                  Y_best     = Y;
                  Rz_best    = Rz;
                  Rqz_best   = Rqz;  
                  perf       = [perf, best];
                  plot(perf, 'r*-');
                  drawnow; 
                end
                para.fval    = best;
              end
            end
          end
        end
      end      
       para.p         =  p(ind(1));
       para.w_between =  w_between(ind(2));
       para.w_within  =  w_within(ind(3));   
       para.alpha     =  alpha(ind(4));
       para.dim       =  dim(ind(5));
       
    case 'ga'
    
      fixed_para{1}=[]; 
      Ipara{1}=[];
      Iparaname{1} =[];
      Cpara{1}=[];
      Cparaname{1} =[];
      
      if ~isempty(p1);
        fixed_para{end+1} = 'p1';
        fixed_para{end+1} =  p1;
      end
      
      if ~isempty(p2);
        fixed_para{end+1} = 'p2';
        fixed_para{end+1} =  p2;
      end
      
      if ~isempty(w_within1);
        fixed_para{end+1} = 'w_within1';
        fixed_para{end+1} =  w_within1;
      end
      
      if ~isempty(w_within2);
        fixed_para{end+1} = 'w_within2';
        fixed_para{end+1} =  w_within2;
      end
      
      if isempty(p1)* isempty(p2)
        if length(p) ==2
          Ipara{end+1} = p; 
          Ipara{end+1} = p; 
          Iparaname{end+1} = 'p1'; 
          Iparaname{end+1} = 'p2'; 
        else
          fixed_para{end+1} = 'p1';
          fixed_para{end+1} =  p;
          fixed_para{end+1} = 'p2';
          fixed_para{end+1} =  p;
        end   
      end
      
      if length(dim)==2;
        Ipara{end+1} = dim;
        Iparaname{end+1} = 'dim'; 
      else
        fixed_para{end+1} = 'dim';
        fixed_para{end+1} =  dim;
      end
      Ipara(1)=[];
      Iparaname(1)=[];
      
      if length(w_between) ==2
        Cpara{end+1} = w_between; 
        Cparaname{end+1} = 'w_between';
      else
        fixed_para{end+1} = 'w_between';
        fixed_para{end+1} =  w_between;
      end
      
      if isempty(w_within1)* isempty(w_within2)
        if length(w_within) ==2
          Cpara{end+1} = w_within;
          Cpara{end+1} = w_within;
          Cparaname{end+1} = 'w_within1';
          Cparaname{end+1} = 'w_within2';
        else
           fixed_para{end+1} = 'w_within1';
           fixed_para{end+1} =  w_within;
           fixed_para{end+1} = 'w_within2';
           fixed_para{end+1} =  w_within;
        end
      end
      
      
      if length(alpha) ==2
        Cpara{end+1} = alpha; 
        Cparaname{end+1} = 'alpha';
      else
        fixed_para{end+1} = 'alpha';
        fixed_para{end+1} =  alpha;
      end
      
      
      Cpara(1)=[];
      Cparaname(1)=[];
      fixed_para(1)=[];
      
      IparaIni=[];
      for jj = 1:length(Iparaname)
        for ii = 1 : 2 : length(ini)-1
          if (strcmp(lower(Iparaname{jj}), lower(ini{ii}))) | (strcmp(lower(Iparaname{jj}), 'p1')&strcmp('p', lower(ini{ii})))| (strcmp(lower(Iparaname{jj}), 'p2')&strcmp('p', lower(ini{ii})))
            IparaIni(jj) = ini{ii+1};
          end
        end
      end
      CparaIni=[];
      for jj = 1:length(Cparaname)
        for ii = 1 : 2 : length(ini)-1
          if (strcmp(lower(Cparaname{jj}), lower(ini{ii}))) | (strcmp(lower(Cparaname{jj}), 'w_within1')&strcmp('w_within', lower(ini{ii})))| (strcmp(lower(Cparaname{jj}), 'w_within2')&strcmp('w_within', lower(ini{ii})))
            CparaIni(jj) = ini{ii+1};
          end
        end
      end
      
      x0 =[IparaIni, CparaIni];
      
      o.opt        = opt;
      o.obj_option = obj_option ;
      o.obj_para   = obj_para;
      o.R          =R;
      
      
      if length(Cpara)*length(Ipara)>0
        I = {o, [Iparaname, Cparaname], fixed_para};
        setting = {'Function2Optimize',     @ACASEIGfun,...
                   'FixedArgument',         I,...
                   'IntegerVariable',       Ipara,...
                   'IntegerType',           'parameter',...
                   'ContinuousVariable',    Cpara,...
                   'PopulationSize',        30,...
                   'MaxGen',                300,...
                   'toel',                  1e-10,...
                   };
        if ~isempty(x0)
          setting{end+1} = 'IniPopulation';
          setting{end+1} =x0;
        end
        GAPobj = GAParaOpt(setting);
        GAPobj = optimize(GAPobj);  
        result = GAPobj.output;
        x = [result.Optimal_Ipara, result.Optimal_Cpara];
      elseif length(Cpara)>0
        I = {o, Cparaname, fixed_para};
        setting = {'Function2Optimize',     @ACASEIGfun,...
                   'FixedArgument',         I,...
                   'ContinuousVariable',    Cpara,...
                   'PopulationSize',        30,... 
                   'MaxGen',                300,...
                   'toel',                  1e-10,...
                   };
        if ~isempty(x0)
          setting{end+1} = 'IniPopulation';
          setting{end+1} =x0;
        end
        GAPobj = GAParaOpt(setting);
        GAPobj = optimize(GAPobj);  
        result = GAPobj.output;
        x = result.Optimal_Cpara;
      elseif length(Ipara)>0
        I =  {o, Iparaname, fixed_para};
        setting = {'Function2Optimize',     @ACASEIGfun,...
                   'FixedArgument',         I,...
                   'IntegerVariable',       Ipara,...
                   'IntegerType',           'parameter',...
                   'PopulationSize',        30,...
                   'MaxGen',                300,...
                   'toel',                  1e-10,...
                   };
        if ~isempty(x0)
          setting{end+1} = 'IniPopulation';
          setting{end+1} =x0;
        end
        GAPobj = GAParaOpt(setting);
        GAPobj = optimize(GAPobj);  
        result = GAPobj.output;
        x = result.Optimal_Ipara;
      else
        x=[];
        I = {o, [], fixed_para};
      end
      [y, X_best, Y_best, para, Rqz_best, Rq] = ACASEIGfun(x, I);
      para.fval = y;
  end
  time = cputime -time;
  
  % Display
  o      = dma(1- Rqz_best);
  o.algo = 'covat2';
  o      = o.run();
  A1      =  Rq(o.perm.row, o.perm.col); 
  A2       = Rqz_best(o.perm.row, o.perm.col);figure()
  subplot(2,3,1)
  imagesc(A1); title('Original Relation Matrix (permute recovered)')
  subplot(2,3,2)
  imagesc(A2); title('Recovered Relation Matrix')
  subplot(2,3,3)
  plot( X_best(:,1), X_best(:,2), 'ro',   Y_best(:,1), Y_best(:,2), 'b*')
  
  o      = dma(1- Rq);
  o.algo = 'covat2';
  o      = o.run();

  A1      =  Rq(o.perm.row, o.perm.col); 
  A2       = Rqz_best(o.perm.row, o.perm.col); 
  subplot(2,3,4)
  imagesc(A1); title('Original Relation Matrix (permute original)')
  subplot(2,3,5)
  imagesc(A2); title('Recovered Relation Matrix')

    
      
end


