function obj = training( obj )
%%
% TRAIN - CoEmbeddings model training
%
%  Obj = train(Obj)
%
%           obj     -   co-embeddings object
%
%  Return:
%           Obj     -   co-embeddings object
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

% -------------------------------------------------------------------------
switch (lower(obj.method))             % Selecting the algorithm model
    case 'ledsvd'
        [obj.X, obj.Y, obj.cputime, obj.model, obj.best, obj.exitflag, obj.svds, obj.Rqz, obj.Rq]...
            = LEDsvd(obj.R, obj.dim, obj.optimisation_option, obj.para_range);
    case 'ledfsvd'
        [obj.X, obj.Y, obj.cputime, obj.model, obj.best, obj.exitflag, obj.svds, obj.Rqz, obj.Rq]...
            = LEDFsvd(obj.R, obj.dim, obj.optimisation_option, obj.para_range);
    case 'ledsvd_beta'
        [obj.X, obj.Y, obj.cputime, obj.model, obj.best, obj.exitflag, obj.svds, obj.Rqz, obj.Rq]...
            = LEDsvd_beta(obj.R, obj.dim, obj.optimisation_option, obj.para_range);
    case 'ledfsvd_beta'
        [obj.X, obj.Y, obj.cputime, obj.model, obj.best, obj.exitflag, obj.svds, obj.Rqz, obj.Rq]...
            = LEDFsvd_beta(obj.R, obj.dim, obj.optimisation_option, obj.para_range);
    case 'ledeig'
        [obj.X, obj.Y, obj.cputime, obj.model, obj.best, obj.exitflag, obj.svds, obj.Rqz, obj.Rq]...
            = LEDeig(obj.R, obj.dim, obj.optimisation_option, obj.para_range, obj.display);
    case 'code'
        
        if isempty(obj.dim)
            obj.dim=2;
        end
        
        obj.cputime = cputime;
        
        if sum(obj.R(:))==1
            [obj.X, obj.Y, out] = CODE(obj.R, obj.dim, obj.CODEoptions);
        else
            [obj.X, obj.Y, out] = CODE(obj.R/sum(obj.R(:)),obj.dim, obj.CODEoptions);
        end
        
        obj.likelihood = out.lik;
        
        obj.cputime = cputime - obj.cputime;
        Rz     = CalRelationXY(obj.X, obj.Y, 'euclidean_sim');
        Rz     =  ACAS_R(Rz, 'quant_int', {'bin', 10});
        R      =  ACAS_R(obj.R, 'quant_int', {'bin', 10});
        o      = dma(1- Rz);
        o.algo = 'covat2';
        o      = o.run();
        A1      =  R(o.perm.row, o.perm.col);
        A2       = Rz(o.perm.row, o.perm.col);
        obj.Rqp =A1;
        obj.Rqpz=A2;
        
        if obj.display
            figure('Visible','Off');
            
            subplot(2,3,1)
            imagesc(A1); title('Original Relation Matrix (permute recovered)')
            subplot(2,3,2)
            imagesc(A2); title('Recovered Relation Matrix')
            subplot(2,3,3)
            plot( obj.X(:,1), obj.X(:,2), 'ro',   obj.Y(:,1), obj.Y(:,2), 'b*')
            o      = dma(1- obj.R);
            o.algo = 'covat2';
            o      = o.run();
            A1      =  obj.R(o.perm.row, o.perm.col);
            A2       = Rz(o.perm.row, o.perm.col);
            subplot(2,3,4)
            imagesc(A1); title('Original Relation Matrix (permute original)')
            subplot(2,3,5)
            imagesc(A2); title('Recovered Relation Matrix')
        end
        
        
    case {'acas', 'acas_svd', 'acassvd'}
        [obj.X, obj.Y, obj.cputime, obj.model, obj.Rqp, obj.Rqpz,obj.best] = ACAS(obj.R, obj.dim, obj.para_range, obj.optimisation_option, obj.ini, obj.display);
        
    case {'acas_eig', 'acaseig'}
        [obj.X, obj.Y, obj.cputime, obj.model] = ACASeig(obj.R, obj.dim, obj.para_range, obj.optimisation_option, obj.ini, obj.display);
    case 'lsi'
        obj.optimisation_range = {'p',0, 'alpha', 0, 'beta', 0};
        obj.optimisation_option = {'optimization', 'grid'};
        [obj.X, obj.Y, obj.cputime, obj.model, obj.Rqp, obj.Rqpz] = ACAS(obj.R, obj.dim, obj.para_range, obj.optimisation_option, obj.ini, obj.display);
    case 'ca'
        obj.para_range = {'p',1, 'alpha', 0.5, 'beta', 1};
        obj.optimisation_option = {'optimization', 'grid'};
        [obj.X, obj.Y, obj.cputime, obj.model, obj.Rqp, obj.Rqpz] = ACAS(obj.R, obj.dim, obj.para_range, obj.optimisation_option, obj.ini, obj.display);
    case 'bgp'
        obj.para_range = {'p',1, 'alpha', 0.5, 'beta', 0};
        obj.optimisation_option = {'optimization', 'grid'};
        [obj.X, obj.Y, obj.cputime, obj.model, obj.Rqp, obj.Rqpz] = ACAS(obj.R, obj.dim, obj.para_range, obj.optimisation_option, obj.ini, obj.display);
    otherwise
        error(['unreliable methods: ' obj.method]);
end

end

