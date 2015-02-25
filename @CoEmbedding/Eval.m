function obj = Eval( obj )
%%
% Eval - CoEmbeddings model objective function evaluation
%
%  Obj = Eval(Obj)
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
switch (obj.method)             % Selecting the algorithm model
    case 'LEDsvd'
         LEDsvd_eval(obj.R, obj.dim, obj.optimisation_option, obj.para_range);
end

end

