function obj = subsasgn(obj, index, val)

  switch index.type
    case '()'
      switch index.subs{:}
        otherwise error('index out of range')
      end
    case '.'
      switch index.subs 
        case 'MyModel'
          obj.MyModel = val;
        case 'FixedArgument'
          obj.FixedArgument =val;
        case 'Ipara'
          obj.Ipara   = val;
        case 'Itype'
          obj.Itype   = val;
        case 'Cpara'
          obj.Cpara   = val;
        case 'PopulationSize'
          obj.PopulationSize = val;
        case 'CrossoverRate'
          obj.CrossoverRate = val;
        case 'MutationUniformRatio'
          obj.MutationUniformRatio  = val; 
        case 'IniPopulation'
          obj.IniPopulation = val;
        case 'MaxGen'
          obj.MaxGen        = val;
        case 'toel'
          obj.toel =val;
        case 'output'
          obj.output   =val;
        case 'display'
          obj.display               =  val;
        otherwise
          error('invalid field name')
      end
    otherwise
      error('undefined assignment method')
  end

end
