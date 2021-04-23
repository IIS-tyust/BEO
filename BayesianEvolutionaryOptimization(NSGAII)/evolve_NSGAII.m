function  [realIndividual,Population,fitnessP,EIP,FrontNo,CrowdDis] = evolve_NSGAII(dmodel,Archive_train,Population,fitnessP,EIP,Boundary,FrontNo,CrowdDis)
w = 0;
wmax = 20;
[N,D] = size(Population);
m = size(Archive_train,2) - D;
 f_min = min(Archive_train(:,D+1:D+m),[],1);
%  Coding = 'Real';
    %% Optimization
    fitnessO = zeros(N,m);
     ER = fitnessO; s = ER; EIO = fitnessO;
    while w <= wmax
        MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
        Offspring  = GA(Population(MatingPool',:),Boundary);
%         Offspring = P_generator(Population(MatingPool',:),Boundary,Coding,N);

        for i = 1: size(Offspring,1)
            for t= 1:m
                [fitnessO(i,t),~,ER(i,t)] = predictor(Offspring(i,:),dmodel(t));
                s(i,t) = sqrt(ER(i,t));
                EIO(i,t) = (f_min(1,t) - fitnessO(i,t)) * Gaussian_CDF((f_min(1,t)-fitnessO(i,t))./ s(i,t))...
                    + s(i,t) * Gaussian_PDF((f_min(1,t)-fitnessO(i,t))./ s(i,t));
            end
        end
        fitness = [fitnessP;fitnessO];
        Population = [Population; Offspring];
        EI = [EIP; EIO];
%%    环境选择 fitness and EI space
        [Population,fitnessP,Index,FrontNo,CrowdDis] = EnvironmentalSelection(Population,fitness,N);
          EIP = EI(Index,:);
          
%         [Population,EIP,~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,EI,N);
%          size(Population(find(FrontNo==1)),1)
        w = w + 1;
    end
    %%  first front(环境选择出来的第一层个体再按EI准则选择个体)
      PopulationFirst = Population(find(FrontNo==1),:); 
      EIfirst = EIP(find(FrontNo==1),:);
     [FrontNum,~] = NDSort(EIfirst,inf);
     realIndividual = PopulationFirst(find(FrontNum==1),:); 
     realIndividual = unique(realIndividual,'rows');
    %%  全程EI空间，直接选出第一层个体
%      realIndividual  = Population(find(FrontNo==1),:); 
%      realIndividual = unique(realIndividual,'rows');
     %% m+1
%      PopulationFirst = Population(find(FrontNo==1),:);
%      EIfirst = EIP(find(FrontNo==1),:);
%      [FrontNum,~] = NDSort(EIfirst,inf); 
%      PopulationEI = PopulationFirst(find(FrontNum==1),:); 
%      EIfirstf = EIfirst(find(FrontNum==1),:);
%      EInum = size(PopulationEI,1)
%      if EInum > m+1
%          EI_mid = zeros(1,m);
%          for k = 1:m
%              [~,ind] = min(EIfirstf(:,k));
%              EI_mid(1,k) = (1/2)*max(EIfirstf(:,k));
% %              EI_mid(1,k) = (1/2)*(max(EIfirstf(:,k))-min(EIfirstf(:,k)));
%              realIndividual(k,:) = PopulationEI(ind,:);
%          end
%          [~,index] = min(dist(EIfirstf,EI_mid'));
%          realIndividual(m+1,:) = PopulationEI(index,:);
%      else
%          realIndividual = PopulationEI;
%      end
%      realIndividual = unique(realIndividual,'rows');

     %%  全程EI空间 从第一层选择m+1个个体
%      PopulationFirst = Population(find(FrontNo==1),:);
%      EIfirstf = EIP(find(FrontNo==1),:);
%      firstnum = size(PopulationFirst,1)
%      if firstnum > m+1
%          for k = 1:m
%              [~,ind] = min(EIfirstf(:,k));
%              EI_mid(1,k) = (1/2)*max(EIfirstf(:,k));
%            %  EI_mid(1,k) = (1/2)*(max(EIfirstf(:,k))-min(EIfirstf(:,k)));
%              realIndividual(k,:) = PopulationFirst(ind,:);
%          end
%          [~,index] = min(dist(EIfirstf,EI_mid'));
%          realIndividual(m+1,:) = PopulationFirst(index,:);
%      else
%          realIndividual = PopulationFirst;
%      end
%      realIndividual = unique(realIndividual,'rows');
end