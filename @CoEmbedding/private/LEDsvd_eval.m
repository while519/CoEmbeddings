function [  ] = LEDsvd_eval( R, dim, options, para_range )
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

activate_variable = [];
% --------------------assign parameter values------------------------------
for ii = 1 : 2 : length(para_range)-1
    switch lower(para_range{ii})
        case 'eta1'
            if isscalar(para_range{ii+1})
                eta1 = para_range{ii+1};
            else
                activate_variable = [activate_variable 1];
                eta1_lower = para_range{ii+1}(1);
                eta1_upper = para_range{ii+1}(2);
            end
        case 'eta2'
            if isscalar(para_range{ii+1})
                eta2 = para_range{ii+1};
            else
                activate_variable = [activate_variable 2];
                eta2_lower = para_range{ii+1}(1);
                eta2_upper = para_range{ii+1}(2);
            end
        case 'alpha'
            if isscalar(para_range{ii+1})
                alpha = para_range{ii+1};
            else
                activate_variable = [activate_variable 3];
                alpha_lower = para_range{ii+1}(1);
                alpha_upper = para_range{ii+1}(2);
            end
        case 'beta'
            if isscalar(para_range{ii+1})
                beta = para_range{ii+1};
            else
                activate_variable = [activate_variable 4];
                beta_lower = para_range{ii+1}(1);
                beta_upper = para_range{ii+1}(2);
            end
        case 'zeta'
            if isscalar(para_range{ii+1})
                zeta = para_range{ii+1}
            else
                activate_variable = [activate_variable 5];
                zeta_lower = para_range{ii+1}(1);
                zeta_upper = para_range{ii+1}(2);
            end
        case 'interval'
            interval = para_range{ii+1};
        otherwise
            error( [ 'unknown model variables: ' para_range{ii} ] )
    end
end

% ---------------------optimisation options--------------------------------
for ii = 1 : 2 : length(options)-1
    switch lower(options{ii})
        case {'evaluation', 'eval'}
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
if ~(exist('zeta_lower','var') || exist('zeta', 'var'))
    switch lower(opt)
        case '2d'
            if ~isscalar(activate_variable)
                error('More than on 2 dimensional plot');
            else
                switch activate_variable
                    case 1
                        S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable eta1']);
                        disp(S);
                        x = eta1_lower : interval: eta1_upper;
                        y = LEDsvd_obj([x', repmat([eta2, alpha, beta],length(x),1)], R, dim, obj_type, obj_para);
                        figure;
                        plot(x,y,'.');
                    case 2
                        S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable eta2']);
                        disp(S);
                        x = eta2_lower : interval: eta2_upper;
                        y = LEDsvd_obj([repmat(eta1, length(x),1), x', repmat([ alpha, beta],length(x),1)], R, dim, obj_type, obj_para);
                        figure;
                        plot(x,y,'.');
                    case 3
                        S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable alpha']);
                        disp(S);
                        x = alpha_lower : interval: alpha_upper;
                        y = LEDsvd_obj([repmat([eta1,eta2],length(x),1), x', repmat( beta,length(x),1)], R, dim, obj_type, obj_para);
                        
                        figure
                        plot(x,y,'.');
                    case 4
                        S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable beta']);
                        disp(S);
                        x = beta_lower : interval: beta_upper;
                        y = LEDsvd_obj([ repmat([eta1, eta2, alpha],length(x),1), x'], R, dim, obj_type, obj_para);
                        figure
                        plot(x,y,'.');
                    otherwise
                        error( 'unknow error in the switch case in LEDsvd_eval');
                end
            end
        case '3d'
            if length(activate_variable) > 2
                error('more than 2 variable has been activate');
                
            elseif ismember(activate_variable, [1,2])           % eta1 and eta2
                S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable eta1, eta2']);
                        disp(S);
                x = eta1_lower: interval : eta2_upper;
                y = eta2_lower: interval : eta2_upper;
                [X, Y] = meshgrid(x,y);
                Z = LEDsvd_obj([ X(:),Y(:),repmat( [alpha beta],length(X(:)),1)], R, dim, obj_type, obj_para);
                
                Z = reshape(Z,size(X));
                figure
                mesh(X,Y,Z);
                figure
                contour(X,Y,Z);
                
            elseif ismember([3,4],activate_variable)            % alpha and beta
                S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variable alpha, beta']);
                        disp(S);
                x = alpha_lower: interval : alpha_upper;
                y = beta_lower: interval : beta_upper;
                [X, Y] = meshgrid(x,y);
                Z = LEDsvd_obj([ repmat( [eta1 eta2],length(X(:)),1), X(:),Y(:)], R, dim, obj_type, obj_para);
                
                Z = reshape(Z,size(X));
                figure
                mesh(X,Y,Z);
                figure
                contour(X,Y,Z);
                
            elseif activate_variable == 3                        % double alpha scaling
                S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variables muilti_alpha']);
                        disp(S);
                x = alpha_lower: interval : alpha_upper;
                [X, Y] = meshgrid(x,x);
                Z = LEDFsvd_obj([ repmat( [eta1 eta2],length(X(:)),1), X(:),Y(:),repmat( beta, length(X(:)),1)], R, dim, obj_type, obj_para);
                
                Z = reshape(Z,size(X));
                figure
                mesh(X,Y,Z);
                figure
                contour(X,Y,Z);
            elseif activate_variable == 4                       % double beta scaling
                                S=sprintf(['Plot for LEDsvd model with ' obj_type ' for variables muilti_beta']);
                        disp(S);
                x = beta_lower: interval : beta_upper;
                [X, Y] = meshgrid(x,x);
                Z = LEDsvd_obj([ repmat( [eta1 eta2, alpha],length(X(:)),1), X(:),Y(:)], R, dim, obj_type, obj_para);
                
                Z = reshape(Z,size(X));
                figure
                mesh(X,Y,Z);
                figure
                contour(X,Y,Z);
                
            end
        otherwise
            error( ['unrecognised optimisation method ' opt ])
    end
else
    %~~~~~~~~~~~~~~~~~~~~~~~~~~zeta parameter model~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    switch lower(opt)
        case '2d'
            if ~isscalar(activate_variable)
                error('More than on 2 dimensional plot');
            else
                switch activate_variable
                    case 1
                        S=sprintf(['Plot for LEDsvdz model with ' obj_type ' for variable eta1']);
                        disp(S);
                        x = eta1_lower : interval: eta1_upper;
                        y = LEDsvdz_obj([x', repmat([eta2, zeta, beta],length(x),1)], R, dim, obj_type, obj_para);
                        figure
                        plot(x,y,'.');
                    case 2
                        S=sprintf(['Plot for LEDsvdz model with ' obj_type ' for variable eta2']);
                        disp(S);
                        x = eta2_lower : interval: eta2_upper;
                        y = LEDsvdz_obj([repmat(eta1, length(x),1), x', repmat([ zeta, beta],length(x),1)], R, dim, obj_type, obj_para);
                        figure
                        plot(x,y,'.');
                    case 5
                        S=sprintf(['Plot for LEDsvdz model with ' obj_type ' for variable zeta']);
                        disp(S);
                        x = zeta_lower : interval: zeta_upper;
                        y = LEDsvdz_obj([repmat([eta1,eta2],length(x),1), x', repmat( beta,length(x),1)], R, dim, obj_type, obj_para);
                        figure
                        plot(x,y,'.');
                    case 4
                        S=sprintf(['Plot for LEDsvdz model with ' obj_type ' for variable beta']);
                        disp(S);
                        x = beta_lower : interval: beta_upper;
                        y = LEDsvdz_obj([ repmat([eta1, eta2, zeta],length(x),1), x'], R, dim, obj_type, obj_para);
                        figure
                        plot(x,y,'.');
                    otherwise
                        error( 'unknow error in the switch case in LEDsvd_eval');
                end
            end
    end
    
end

