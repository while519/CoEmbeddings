function Population = CreateIC(GenomeLength,FitnessFcn,options)

%CreateIC Creates the initial population for genetic algorithm, including a mixture of integer and continuous variables.
% For continuous variables, the MATLAB function "gacreationlinearfeasible" is used.
% For integer variables, random generatation is used.

%   Example:
%     options = gaoptimset('CreationFcn',@CreateIC);

%   T. Mu June 2011

  load('GAinfor')
  GenomeLengthI = length(Ipara);
  GenomeLengthC = length(Cpara);
  
  totalPopulation = sum(options.PopulationSize);
  initPopProvided = size(options.InitialPopulation,1);
  individualsToCreate = totalPopulation - initPopProvided;

  if GenomeLengthC > 0 && GenomeLengthI >0    
    %Treat both integer and continuous variables as continuous first for generating population.
    optionsC = options; optionsC.PopulationType = 'doubleVector';
    Population = gacreationlinearfeasible(GenomeLength,[],optionsC);
    %Replace generated population for ineger with integer population
    for ii=1:GenomeLengthI
      tt = Ipara{ii}; Iopt = tt(1):tt(2); Nopt = length(Iopt);
      for jj=1:individualsToCreate
        ind = randperm(Nopt);
        pop(jj, ii) = Iopt(ind(1));
      end
    end
    Population(initPopProvided+1:end, 1:GenomeLengthI) = pop;
    %When user-defined initial population is not provided, GA tool automatically generate one individual, which can be continuous for interger variable
    Population(1:initPopProvided,1:GenomeLengthI) = round(Population(1:initPopProvided,1:GenomeLengthI));      
  elseif GenomeLengthC == 0 && GenomeLengthI >0    
    Population = zeros(totalPopulation,GenomeLengthI);
    if initPopProvided > 0  
       Population(1:initPopProvided,:) = round(options.InitialPopulation);
    end
    for ii=1:GenomeLengthI
      tt = Ipara{ii}; Iopt = tt(1):tt(2); Nopt = length(Iopt);
      for jj=1:individualsToCreate
        ind = randperm(Nopt);
        pop(jj, ii) = Iopt(ind(1));
      end
    end
    Population(initPopProvided+1:end,:)=pop;    
  elseif GenomeLengthI==0 && GenomeLengthC >0    
    optionsC = options; optionsC.PopulationType = 'doubleVector';
    Population = gacreationlinearfeasible(GenomeLength,[],optionsC);    
  end  
  
end

