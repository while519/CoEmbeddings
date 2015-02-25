function [ Xq, Yq, time, para, best, exitflag, svds, Rqz, Rq ] = LEDsvd_beta( R, dim, options, para_range )
%%
% LEDFSVD - Runing the LEDF model using svd and alpha scaling
%
%  [ Xq, Yq, time, para, best ] = LEDFsvd( R, dim, options, para_range )
%
%           R           -   (MxN) input relational matrix
%           dim         -   (scalor) tageted embedded dimension
%           options     -   for xxx score setting
%           para_range  -   parameters( alpha, beta, eta1, eta2) range
%
%  Return:
%           Xq,Yq     -   (M,Nxdim) row and column embeddings
%           time      -   running time for the model
%           para      -   (struct) parameters for best matched model
%           best      -   best criterion function( score) value
%
%%
%  Author   : Yu Wu
%            University of Liverpool
%            Electrical Engineering and Electronics
%            Brownlow Hill, Liverpool L69 3GJ
%            yuwu@liv.ac.uk
%  Last Rev : Friday, January 17, 2014 (GMT) 10:03 AM
%  Tested   : Matlab_R2014b
%
%%
% Copyright notice: You are free to modify, extend and distribute
%       this code granted that the author of the original code is
%       mentioned as the original author of the code.

% --------------------assign parameter values------------------------------
for ii = 1 : 2 : length(para_range)-1
    switch lower(para_range{ii})
        case 'eta1'
            if isscalar(para_range{ii+1})
                eta1_upper = para_range{ii+1};
                eta1_lower = para_range{ii+1};
            else
                eta1_lower = para_range{ii+1}(1);
                eta1_upper = para_range{ii+1}(2);
            end
        case 'eta2'
            if isscalar(para_range{ii+1})
                eta2_lower = para_range{ii+1};
                eta2_upper = para_range{ii+1};
            else
                eta2_lower = para_range{ii+1}(1);
                eta2_upper = para_range{ii+1}(2);
            end
        case 'alpha'
            if isscalar(para_range{ii+1})
                alpha_lower = para_range{ii+1};
                alpha_upper = para_range{ii+1};
            else
                alpha_lower = para_range{ii+1}(1);
                alpha_upper = para_range{ii+1}(2);
            end
        case 'beta'
            if isscalar(para_range{ii+1})
                beta_lower = para_range{ii+1};
                beta_upper = para_range{ii+1};
            else
                beta_lower = para_range{ii+1}(1);
                beta_upper = para_range{ii+1}(2);
            end
        case 'zeta'
            if isscalar(para_range{ii+1})
                zeta_lower = para_range{ii+1};
                zeta_upper = para_range{ii+1};
            else
                zeta_lower = para_range{ii+1}(1);
                zeta_upper = para_range{ii+1}(2);
            end
        otherwise
            error( [ 'unknown model variables: ' para_range{ii} ] )
    end
end

% ---------------------optimisation options--------------------------------
for ii = 1 : 2 : length(options)-1
    switch lower(options{ii})
        case {'optimisation', 'optimization'}
            opt = options{ii+1};
        case {'obj_option', 'score'}
            obj_type = options{ii+1};
        case 'obj_para'
            obj_para = options{ii+1};
        otherwise
            error( ['unknown setting for optimisation type ' options{ii}])
    end
end

% _____________________________alpha parameter model_______________________
if ~exist('zeta_lower','var')
    switch lower(opt)
        case 'ga'
            clc
            S=sprintf(['It is current running LEDsvd_beta' '.....eq_alpha']);
            disp(S);
            options1 = gaoptimset('PopulationSize',52,'UseParallel', false,...
                'Vectorized', 'on','Display','iter',...
                'PlotFcns',@gaplotbestf...
                );
            time = cputime;
            [x,fval,exitflag] = ga(@(x) LEDsvd_obj(x, R, dim, obj_type, obj_para),3+dim,[],[],[],[],...
                [eta1_lower, eta2_lower, alpha_lower, repmat(beta_lower,1,dim)],...
                [eta1_upper, eta2_upper, alpha_upper, repmat(beta_upper,1,dim)],...
                [],[],options1);
            para.eta1 = x(1);
            para.eta2 = x(2);
            para.alpha = x(3);
            para.beta = x(4:end);
            time = cputime-time;
            best = fval;
            [~,Xq,Yq,svds, Rqz, Rq] = LEDsvd_model(x, R, dim, obj_type, obj_para);
        otherwise
            error( ['unrecognised optimisation method ' opt ])
    end
else
%~~~~~~~~~~~~~~~~~~~~~~~~~~zeta parameter model~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    switch lower(opt)                           
        case 'ga'
            clc
            S=sprintf(['It is current running LEDsvdz_beta' '.....eq_zeta']);
            disp(S);
            options1 = gaoptimset('PopulationSize',52,'UseParallel', false,...
                'Vectorized', 'on','Display','iter',...
                'PlotFcns',@gaplotbestf...
                );
            time = cputime;
            [x,fval,exitflag] = ga(@(x) LEDsvdz_obj(x, R, dim, obj_type, obj_para),dim+3,[],[],[],[],...
                [eta1_lower, eta2_lower, zeta_lower, repmat(beta_lower,1,dim)],...
                [eta1_upper, eta2_upper, zeta_upper, repmat(beta_upper,1,dim)],...
                [],[],options1);
            para.eta1 = x(1);
            para.eta2 = x(2);
            para.zeta = x(3);
            para.beta = x(4:end);
            time = cputime-time;
            best = fval;
            [~,Xq,Yq,svds, Rqz, Rq] = LEDsvdz_model(x, R, dim, obj_type, obj_para);
        otherwise
            error( ['unrecognised optimisation method ' opt ])
    end
end

end

