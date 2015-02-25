function val = subsref(obj, index)

  switch index.type
    case '()'
      switch index.subs{:}
        otherwise error('index out of range')
      end
    case '.'
      switch index.subs
        case 'MyModel'
          val = obj.MyModel;
        case 'FixedArgument'
          val = obj.FixedArgument;
        case 'Ipara'
          val = obj.Ipara    ;
        case 'Cpara'
          val = obj.Cpara    ;
        case 'Itype'
          val = obj.Itype;
        case 'PopulationSize'
          val = obj.PopulationSize;
        case 'CrossoverRate'
          val = obj.CrossoverRate;
        case 'MutationUniformRatio'
          val = obj.MutationUniformRatio; 
        case 'IniPopulation'
          val = obj.IniPopulation  ;
        case 'MaxGen'
          val = obj.MaxGen         ;
        case 'toel'
          val = obj.toel;
        case 'output'
          val  = obj.output;
        case 'display'
          display =  I{ii+1};
        otherwise
          error('invalid field name')
      end
    otherwise
        error('undefined assignment method')
  end
  
end

%------------------------------------------------