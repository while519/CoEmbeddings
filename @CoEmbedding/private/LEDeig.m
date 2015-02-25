function [ Xq, Yq, time, para, best, exitflag, svds, Rqz, Rq ] = LEDeig( R, dim, options, para_range, display )
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
        otherwise
            error( [ 'unknown model variables for LEDeig: ' para_range{ii} ] )
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

switch lower(opt)
    case 'ga'
        clc
        S=sprintf(['It is current running LEDeig' '.....eq_alpha']);
        disp(S);
        if ~isequal(eta1_upper, eta1_lower) || ~isequal(eta2_upper,eta2_lower) || ...
                ~isequal(alpha_lower,alpha_upper) || ~isequal(beta_lower, beta_upper)
            options1 = gaoptimset('PopulationSize',104,'UseParallel', false,...
                'Vectorized', 'on','Display','iter',...
                'PlotFcns',@gaplotbestf...
                );
            time = cputime;
            [x,fval,exitflag] = ga(@(x) LEDeig_obj(x, R, dim, obj_type, obj_para),4,[],[],[],[],...
                [eta1_lower, eta2_lower, alpha_lower, beta_lower],...
                [eta1_upper, eta2_upper, alpha_upper, beta_upper],...
                [],[],options1);
            time = cputime-time;
            best = fval;
        else
            x = [eta1_lower, eta2_lower, alpha_lower, beta_lower];
            time = 0;
            best = 0;
            exitflag = -100;
        end
        
        para.eta1 = x(1);
        para.eta2 = x(2);
        para.alpha = x(3);
        para.beta = x(4);
        [~,Xq,Yq,svds, Rqz, Rq] = LEDeig_model(x, R, dim, obj_type, obj_para);
        
    case 'grid'
        disp(['It is current running LEDeig' '..... eq_alpha with grid search']);
        figure;
        if ~isequal(eta1_upper, eta1_lower) || ~isequal(eta2_upper,eta2_lower) || ...
                ~isequal(alpha_lower,alpha_upper) || ~isequal(beta_lower, beta_upper)
            Rq = ACAS_R(R, obj_type, obj_para);
            best = 1e10;
            perf = [];
            time = cputime;
            for eta1 = eta1_lower:0.5:eta2_upper           % linespace
                Dr=sum(R,2);             % row sum diagonal matrix     column vector
                Rx=bsxfun(@times,R,Dr.^eta1);
                Pc=sum(Rx,1);
                
                T2 = bsxfun(@rdivide, Rx, Pc);
                for eta2 =  eta2_lower:0.5:eta2_upper
                    Dc=sum(R,1);             % column sum diagonal matrix  row vector
                    Ry=bsxfun(@times,R,Dc.^(eta2));
                    Qr=sum(Ry,2);
                    
                    T1 = bsxfun(@rdivide, Ry, Qr);
                    
                    [U,D]=eigs(T1*T2',dim+1);                   % S include the nonnegative eigen values in decresing order
                    [ U, B] = EigSort( U, D, 'descend');
                    
                    p= 2:dim+1;
                    BB=B(p);                      % the k second largest eigenvalue vector
                    sigma = BB/BB(1);
                    for bta = beta_lower:0.1:beta_upper
                        %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        Xp = bsxfun(@times,(sigma.^bta)',U(:,p));
                        gamma = sqrt(1./(diag(U(:,p)'*bsxfun(@times, Qr, U(:,p)))));
                        X = bsxfun(@times, Xp, gamma');
                        for alpha = alpha_lower:0.1:alpha_upper
                            
                            Yp = bsxfun(@rdivide,Rx'*X*alpha,Pc');
                            Y = bsxfun(@rdivide, Yp, sqrt(BB'));
                            
                            [score, ~, ~] = KNN_obj(Rq, X, Y, obj_type, obj_para);
                            if score < best
                                best       = score;
                                x        = [eta1,eta2,alpha,bta];
                            end
                            
%++++++++++++++++++++++++++++++++++ Display +++++++++++++++++++++++++++++++                         
                            perf       = [perf, best];
                            
                            if display
                                plot(perf, 'r.-');
                                drawnow;
                            end
%--------------------------------------------------------------------------                            
                        end
                    end
                end
            end
            time = cputime - time;
        else
            x = [eta1_lower, eta2_lower, alpha_lower, beta_lower];
            time = 0;
            best = 0;
            exitflag = -100;
        end
        
        para.eta1 = x(1);
        para.eta2 = x(2);
        para.alpha = x(3);
        para.beta = x(4);
        
        [~,Xq,Yq,svds, Rqz, Rq] = LEDeig_model(x, R, dim, obj_type, obj_para);
        
    otherwise
        error( ['unrecognised optimisation method ' opt ])
end
end


