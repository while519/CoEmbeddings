function [X_best, Y_best, time, para,  permuted_Rq, permuted_Rqz,best] = ACAS(R, dim, options, matrix_match, ini, display)
% For ga optimisation, the dim, p, alpha of ACAS is either defined as a fixed number or a
% range. For beta, it is either defined as [] (to be searched from [0,1]) or fixed vector with length
% equal to sim_order.
% For grid search, the dim, p, alpha, beta of ACAS is either defined as a fixed number or a
% set of candidate parameters, because of grid search.
  permuted_Rq=[];
  permuted_Rqz=[];
  
 
  obj_para = [];
  obj_option = 'quant_integer';
  opt = 'ga';
  opt_Rz = 'euclidean_sim';
  
  for ii = 1 : 2 : length(matrix_match)-1
    switch lower(matrix_match{ii})
      case 'estimate_r'
        opt_Rz = matrix_match{ii+1};
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
      p           = [-1, 0, 1, 2, 3];
      alpha       = 0:0.1:5;
      beta        = 0:0.2:2;  
      if isempty(dim)
        temp = min(size(R,1), size(R,2));
        dim = 2: floor((temp-2)/20)+1 : temp ;
      end
    case 'ga'
      p           = [-1, 10];
      p1          = [-1, 10];
      p2          = [-1, 10];
      alpha       = [0, 5];
      beta        = [0, 5];
      if isempty(dim)
        dim = [2, min(size(R,1), size(R,2))];
      end
  end
  
  marker =1;
  for ii = 1 : 2 : length(options)-1
    switch lower(options{ii})
      case 'dim'
        if length(dim)>1
          dim = options{ii+1};
        end
      case 'p'
        p = options{ii+1};
        marker =0;        
      case 'p1'
        p1 = options{ii+1};  
      case 'p2'
        p2 = options{ii+1};        
      case 'alpha'
        alpha  = options{ii+1};
      case 'beta'
        beta = options{ii+1};
      otherwise
        error( [ 'unknown model variables: ' options{ii} ] );
    end
  end
  
  if isempty(ini)
    if marker==0
      ini= {'dim', 2, 'p', 1, 'alpha', 0.5, 'beta', 1};
    else
      ini= {'dim', 2, 'p1', 1, 'p2', 1, 'alpha', 0.5, 'beta', 1};
    end
  end
  
  
  fprintf('ACAS Model Optimisation:\n')
  time = cputime;  
  switch  lower(opt)
    case {'grid_search', 'grid'}   
      p0=-10;
      best = 1e10;
      perf = [];     
      perf_para = [];
      switch obj_option
        case {'quant_integer', 'quant_int', 'quant_adapt',...
              'fast_quant_integer', 'fast_quant_int', 'fast_quant_adapt'}
          Rq      =  ACAS_R(R, obj_option, obj_para);
      end
      for ii=1:length(p)
        for jj=1:length(alpha)
          for kk=1:length(beta)  
            for tt =1:length(dim)
              if p(ii)~=p0
                
                [X, Y, SVDres] = ACASmodel(R,  opt, dim(tt),  p(ii),  alpha(jj),  beta(kk)); 
                
                p0 = p(ii);
              else
                [X, Y] = ACASmodel(R,  opt, dim(tt),  p(ii),  alpha(jj),  beta(kk), SVDres);  
              end
              
              Rz     = CalRelationXY(X, Y, opt_Rz, 'mean'); 
              
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
              switch obj_option
                case {'quant_integer', 'quant_int', 'quant_adapt',...
                    'fast_quant_integer', 'fast_quant_int', 'fast_quant_adapt'}
                  [score, Rqz]  = ACAS_obj(Rq, Rz, obj_option, obj_para);
                otherwise
                  [score, Rqz, Rq]  = ACAS_obj(R, Rz, obj_option, obj_para);
              end
             
              
              if score < best
                best       = score;
                ind        = [ii,jj,kk,tt];
                X_best     = X;
                Y_best     = Y;
                Rz_best    = Rz;
                Rqz_best   = Rqz; 
                best_para.dim = dim(tt);
                best_para.p   = p(ii);
                best_para.alpha    =alpha(jj);
                best_para.beta     = beta(kk);
              end              
              perf       = [perf, best];
              if length(p)*length(alpha)*length(beta)*length(dim)>1
                if display
                  plot(perf, 'r.-');
                  drawnow; 
                end
              end
              perf_para{end+1}= best_para; 
              save 'ACAS_grid_iteration' perf perf_para
              para.fval    = best;              
            end
          end
        end
      end      
      best
       para.p     =  p(ind(1));
       para.alpha =  alpha(ind(2));
       para.beta  =  beta(ind(3));   
       para.dim   =  dim(ind(4));
       
    case 'ga'
    
      fixed_para{1}=[]; 
      Ipara{1}=[];
      Iparaname{1} =[];
      Cpara{1}=[];
      Cparaname{1} =[];
      
                 
      
      if marker==0
        if length(p)==2;
          Ipara{end+1} = p;
          Iparaname{end+1} = 'p'; 
        else
          fixed_para{end+1} = 'p';
          fixed_para{end+1} =  p;
        end
      else
        if length(p1)==2;
          Ipara{end+1} = p1;
          Iparaname{end+1} = 'p1'; 
        else
          fixed_para{end+1} = 'p1';
          fixed_para{end+1} =  p1;
        end      
        if length(p2)==2;
          Ipara{end+1} = p2;
          Iparaname{end+1} = 'p2'; 
        else
          fixed_para{end+1} = 'p2';
          fixed_para{end+1} =  p2;
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
      
      if length(alpha) ==2
        Cpara{end+1} = alpha; 
        Cparaname{end+1} = 'alpha';
      else
        fixed_para{end+1} = 'alpha';
        fixed_para{end+1} =  alpha;
      end
      
      if length(beta) ==2
        Cpara{end+1} = beta;
        Cparaname{end+1} = 'beta';
      else
         fixed_para{end+1} = 'beta';
         fixed_para{end+1} =  beta;
      end

      
      Cpara(1)=[];
      Cparaname(1)=[];
      fixed_para(1)=[];
      
      IparaIni=[];
      for jj = 1:length(Iparaname)
        for ii = 1 : 2 : length(ini)-1
          %if (strcmp(lower(Iparaname{jj}), lower(ini{ii}))) | (strcmp(lower(Iparaname{jj}), 'p1')&strcmp('p', lower(ini{ii})))| (strcmp(lower(Iparaname{jj}), 'p2')&strcmp('p', lower(ini{ii})))
          if strcmp(lower(Iparaname{jj}), lower(ini{ii}))
           IparaIni(jj) = ini{ii+1};
          end
        end
      end
      CparaIni=[];
      for jj = 1:length(Cparaname)
        for ii = 1 : 2 : length(ini)-1
          if strcmp(lower(Cparaname{jj}), lower(ini{ii}))
            CparaIni(jj) = ini{ii+1};
          end
        end
      end
      
      x0 =[IparaIni, CparaIni];
      
      
      o.opt        = opt;
      o.obj_option = obj_option ;
      o.obj_para   = obj_para;
      o.R          = R;
      o.opt_Rz     = opt_Rz;
      
      
      if length(Cpara)*length(Ipara)>0
        I = {o, [Iparaname, Cparaname], fixed_para};
        setting = {'Function2Optimize',     @ACASfun,...
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
        setting = {'Function2Optimize',     @ACASfun,...
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
        setting = {'Function2Optimize',     @ACASfun,...
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
      [y, X_best, Y_best, para, Rqz_best, Rq] = ACASfun(x, I);
      para.fval = y;
  end
  time = cputime -time;
  o      = dma(1- Rqz_best);
  o.algo = 'covat2';
  o      = o.run();
  A1      =  Rq(o.perm.row, o.perm.col); 
  A2       = Rqz_best(o.perm.row, o.perm.col);
  permuted_Rq =A1;
  permuted_Rqz = A2;  
  % Display
  if 	display
    figure('Visible','Off');

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

    
      
end


