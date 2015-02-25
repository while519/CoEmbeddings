function mutationChildren = MutationIC(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate)
%   Mutation for a mixture of integer and continuous variables.
%   For integer variables,  uniform multi-point mutation is used "mutationuniform"
%   For continuous variables, 'mutationadaptfeasible' is used
%   
%
%   Example:
%     options = gaoptimset('MutationFcn', @MutationIC); 
%
%   T. Mu June 2011

  if nargin < 8 || isempty(mutationRate)
    mutationRate = 0.01; % default mutation rate
  end

  load('GAinfor')
  GenomeLengthI = length(Ipara);
  GenomeLengthC = length(Cpara);
  
  if GenomeLengthC > 0 && GenomeLengthI >0 
    thisPopulationI = thisPopulation(:,1:GenomeLengthI);
    thisPopulationC = thisPopulation(:,GenomeLengthI+1:end);
    mutationChildrenI = mutationInteger(parents,options,GenomeLengthI,FitnessFcn,state,thisScore,thisPopulationI,mutationRate);
    
    optionsC = options; optionsI.PopulationType = 'doubleVector';
    optionsC.LinearConstr.lb(1:GenomeLengthI)=[];
    optionsC.LinearConstr.ub(1:GenomeLengthI)=[];
    mutationChildrenC = mutationadaptfeasible(parents,optionsC,GenomeLengthC, FitnessFcn,state,thisScore,thisPopulationC);
    mutationChildren  = [mutationChildrenI, mutationChildrenC];    
  elseif GenomeLengthC == 0 && GenomeLengthI >0
    mutationChildren = mutationInteger(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate);         
  elseif GenomeLengthI == 0 && GenomeLengthC >0
    mutationChildren = mutationadaptfeasible(parents,options,GenomeLength, FitnessFcn,state,thisScore,thisPopulation);    
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
