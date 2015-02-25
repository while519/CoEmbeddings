function obj = CoEmbedding(I)
%%
% COEMBEDDING - Specifying the class of co-embeddings object
%
%  Obj = CoEmbedding(Setting)
%
%           Setting -   user provided cell array
%
%  Return:
%           Obj     -   co-embeddings object
%
% Example:
%
%   o = CoEmbedding(I);
%   I = {'R', R,...
%       'method', 'LEDF',...
%       'dim',2,...
%       'option',option...
%       };
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
if ~exist('I','var')
    error('Wrong input data given to the constructor');
end

%% Default setting
obj.method = 'LEDF';
obj.dim = 2;
obj.para_range = [];
obj.optimisation_option = {'optimisation', 'ga', 'obj_option', 'quant_knn',...
    'obj_para', {'knn', 20}...
    };
obj.cputime = 1;

obj.display = 1;
obj.ini = [];

%Setting for CODE                                        
  obj.CODEoptions                            =[];  

%% Assign values to objects
if ~nargin
    error('No input values available for the default constructor');
elseif isa(I,'CoEmbedding')        % input is an co-embedding objects
    obj = I;
else
    for ii = 1 : 2 : length(I)-1
        switch lower(I{ii})
            case {'r', 'input'}
                obj.R = I{ii+1};
            case 'method'
                obj.method = I{ii+1};
            case {'model_optimisation', 'model_optimization'}
                obj.optimisation_option = I{ii+1};
            case {'para_range', 'parameter_range'}
                obj.para_range = I{ii+1};
            case {'k', 'dim'}
                obj.dim = I{ii+1};
            case 'code_options'
                obj.CODEoptions               = I{ii+1};
            case {'initialisation','initialization'}
                obj.ini                   = I{ii+1};
            case 'display'
                obj.display               = I{ii+1};
            otherwise
                error( ['unknown input argument: ' I{ii}] );
        end
    end
end

%% Result initialisation
obj.X = [];
obj.Y = [];
obj.best = [];
obj.model = [];
obj.exitflag = -10;
obj.svds = [];
obj.likelihood  = []; %only for CODE output;
obj.Rq = [];
obj.Rqz = [];
obj.Rqp = [];
obj.Rqpz =[];

%% Create co-embedding objects
obj = class(obj, 'CoEmbedding');


end
