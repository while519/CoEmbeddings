function [obj] = GAParaOpt(I)

%Parameter optimization of a given model (solving a minimisation problem)
%Use: o = GAParaOpt(Setting);
%     Setting = {'Function2Optimize', @myfun,...
%                'FixedArgument', 10,
%                'IntegerVariable', {[0,  1], [0, 8], [1, 10]},...
%                'ContinuousVariable', {[0.5,10], [-3, 3], [2.1, 5], [9.9,20]},...
%                'IntegerType', 'parameter' (or 'option'),...
%                }
%function y = myfun(x, FixedArgument);
%Output: o.output
%
%
%-----------------------------------
% T. Mu June 2011
%-----------------------------------
  if ~exist('I', 'var')
    error('wrong input data given to the constuctor')
  end
  %Default Input
  obj.Ipara               = {};
  obj.Itype               = 'parameter';
  obj.Cpara               = {}; 
  obj.MyModel             = [];
  obj.FixedArgument       = [];
  
  %Default GA parameters
  obj.PopulationSize       = 20;
  obj.CrossoverRate        = 0.8;
  obj.MutationUniformRatio = 0.1; %set the ratio for uniform mutation
  obj.IniPopulation        = [];
  obj.MaxGen               = 100;
  obj.toel                 = 1e-10;
  obj.display              = 1;
    
    
  if ~nargin
    error('no default constructor available')
  elseif isa(I, 'GAParaOpt')
    obj = I;
  else
    for ii = 1 : 2 : length(I)-1
      switch I{ii} 
        case 'Function2Optimize'
          obj.MyModel               = I{ii+1};
        case 'FixedArgument'
          obj.FixedArgument         = I{ii+1};
        case 'IntegerVariable'
          obj.Ipara                 = I{ii+1};
        case 'ContinuousVariable'
          obj.Cpara                 = I{ii+1};
        case 'IntegerType'
          obj.Itype                 = I{ii+1};
        case 'PopulationSize'
          obj.PopulationSize        = I{ii+1};
        case 'MaxGen'
          obj.MaxGen                = I{ii+1};
        case 'toel'
          obj.toel                  = I{ii+1};
        case 'CrossoverRate'
          obj.CrossoverRate         = I{ii+1};
        case 'MutationUniformRatio'
          obj.MutationUniformRatio  = I{ii+1};
        case 'IniPopulation'
          obj.IniPopulation         =  I{ii+1};
        case 'display'
          obj.display               =  I{ii+1};
        otherwise
          error( [ 'unknown input argument: ' I{ii} ] )
      end
    end
    %Output
    obj.output              = struct('Optimal_Ipara', [],'Optimal_Cpara', [], 'Optimal_fval',[]);
    obj                     = class(obj,'GAParaOpt'); 
  end
end