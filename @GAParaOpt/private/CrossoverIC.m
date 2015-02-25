function xoverKids = CrossoverIC(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation,ratio)
%   Crossover for a mixture of integer and continuous variables.
%   For integer variables,  'crossoverscattered' is used 
%   For continuous variables, 'crossoverheuristic' is used
%   
%
%   Example:
%     options = gaoptimset('CrossoverFcn', @CrossoverIC); 
%
%   T. Mu June 2011

 if nargin < 7 || isempty(ratio)
    ratio = 1.2;
  end

  load('GAinfor')
  GenomeLengthI = length(Ipara);
  GenomeLengthC = length(Cpara);
  
  if GenomeLengthC > 0 && GenomeLengthI >0 
    thisPopulationI = thisPopulation(:,1:GenomeLengthI);
    thisPopulationC = thisPopulation(:,GenomeLengthI+1:end);
    optionsI = options; 
    optionsI.LinearConstr.lb(GenomeLengthI+1:end)=[];
    optionsI.LinearConstr.ub(GenomeLengthI+1:end)=[];
    switch Itype
      case 'option'
        xoverKidsI  = crossoverscattered(parents,optionsI,GenomeLengthI,FitnessFcn,thisScore,thisPopulationI);
      case 'parameter'
        xoverKidsI  = crossoverheuristic(parents,optionsI,GenomeLengthI,FitnessFcn,thisScore,thisPopulationI);
        xoverKidsI  = round(xoverKidsI);
    end
    optionsC = options;
    optionsC.LinearConstr.lb(1:GenomeLengthI)=[];
    optionsC.LinearConstr.ub(1:GenomeLengthI)=[];
    xoverKidsC  = crossoverheuristic(parents,optionsC,GenomeLengthC,FitnessFcn,thisScore,thisPopulationC,ratio);
    xoverKids  = [xoverKidsI, xoverKidsC];
    
  elseif GenomeLengthC == 0 && GenomeLengthI >0
    switch Itype
      case 'option'
        xoverKids  = crossoverscattered(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation);
      case 'parameter'
        xoverKids  = crossoverheuristic(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation);   
        xoverKids  = round(xoverKids);
    end
  elseif GenomeLengthI == 0 && GenomeLengthC >0
    %xoverKids  = crossoverheuristic(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation,ratio);
     xoverKids  = crossoverintermediate(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation);    
  end  
  
end

%========
function mutationChildren = mutationInteger(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate)


  if nargin < 8 || isempty(mutationRate)
      mutationRate = 0.01; % default mutation rate
  end
  load('GAinfor')
  
  mutationChildren = zeros(length(parents),GenomeLength);
  for i=1:length(parents)
      child = thisPopulation(parents(i),:);
      % Each element of the genome has mutationRate chance of being mutated.
      mutationPoints = find(rand(1,length(child)) < mutationRate);
      for ii=1:length(mutationPoints)
        tt = Ipara{mutationPoints(ii)}; candidate = tt(1):tt(2); 
        candidate( find( candidate == child(mutationPoints(ii)) ) )=[];
        ind = randperm(length(candidate)); 
        child(mutationPoints(ii))= candidate(ind(1));
      end
      mutationChildren(i,:) = child;
  end
  
end
