function [obj] = optimize(obj)

  if isempty(obj.MyModel)
    error('Set the function to be optimized by GA.')
  end
    
  if isempty(obj.Cpara)&&isempty(obj.Ipara)
    error('Set the continuous and (or) integer parameters to be optimized by GA.')
  end

  time = cputime;
  obj = GAopt(obj); 
  fprintf('\n Optimization time: (%2.2fsecs) \n', cputime-time);
 
end