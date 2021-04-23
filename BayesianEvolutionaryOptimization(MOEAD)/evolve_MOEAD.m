function  [realIndividual,Population,fitnessP,EIP,Z] = evolve_MOEAD(dmodel,Archive_train,Population,fitnessP,EIP,Boundary,W,B,T,Z)
    w = 0;
    wmax = 20;
    [N,D] = size(Population);
    m = size(Archive_train,2) - D;
    f_min = min(Archive_train(:,D+1:D+m),[],1);
   %% Parameter setting
    type = 2;
    %%  另外测试
%     Pf = [Population fitnessP EIP];
%     [frontno,~] = NDSort(Pf(:,D+1:D+m),inf);
%     Archive = Pf(find(frontno == 1),:);
%     
    %% Optimization
    Offspring = zeros(N,D); fitnessO = zeros(N,m); ER = fitnessO; s = ER; EIO = fitnessO;
    while w <= wmax
        % For each solution
        for i = 1 : N      
            % Choose the parents
            P = B(i,randperm(size(B,2))); %第i个领域中的参考向量号随机排序

            % Generate an offspring
            Offspring(i,:) = GAhalf(Population(P(1:2),:),Boundary); 
            
%             for i = 1: size(Offspring,1)
                for t= 1:m
                    [fitnessO(i,t),~,ER(i,t)] = predictor(Offspring(i,:),dmodel(t));
                    s(i,t) = sqrt(ER(i,t));
                    EIO(i,t) = (f_min(1,t) - fitnessO(i,t)) * Gaussian_CDF((f_min(1,t)-fitnessO(i,t))./ s(i,t))...
                        + s(i,t) * Gaussian_PDF((f_min(1,t)-fitnessO(i,t))./ s(i,t));
                end
%             end
            % Update the ideal point
            Z = min(Z,fitnessO(i,:));
%             Z = min(Z,EIO(i,:));
            
            % Update the neighbours(目标空间)
            switch type
                case 1
                    % PBI approach
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((fitnessP(P,:)-repmat(Z,T,1)).^2,2));
                    normO   = sqrt(sum((fitnessO(i,:)-Z).^2,2));
                    CosineP = sum((fitnessP(P,:)-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(fitnessO(i,:)-Z,T,1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                case 2
                    % Tchebycheff approach
                    g_old = max(abs(fitnessP(P,:)-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(fitnessO(i,:)-Z),T,1).*W(P,:),[],2);
%                     % Tchebycheff approach(EI空间)
%                     g_old = max(abs(EIP(P,:)-repmat(Z,T,1)).*W(P,:),[],2);
%                     g_new = max(repmat(abs(EIO(i,:)-Z),T,1).*W(P,:),[],2);
                case 3
                    % Tchebycheff approach with normalization
                    Zmax  = max(fitnessP,[],1);
                    g_old = max(abs(fitnessP(P,:)-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                    g_new = max(repmat(abs(fitnessO(i,:)-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                case 4
                    % Modified Tchebycheff approach
                    g_old = max(abs(fitnessP(P,:)-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(fitnessO(i,:)-Z),T,1)./W(P,:),[],2);
            end
            Population(P(g_old>=g_new),:) = repmat(Offspring(i,:),length(P(g_old>=g_new)),1);
            fitnessP(P(g_old>=g_new),:) = repmat(fitnessO(i,:),length(P(g_old>=g_new)),1);
            EIP(P(g_old>=g_new),:) = repmat(EIO(i,:),length(P(g_old>=g_new)),1);
        end
   %%  另外测试
%         Archive = [Archive ;Population fitnessP EIP];
%         [FrontN,~] = NDSort(Archive(:,D+1:D+m),inf);
%         Archive = Archive(find(FrontN==1),:);
        
        w = w + 1;
    end
    
      %%  Archive EI first front
%     [FrontNo,~] = NDSort(Archive(:,D+m+1:D+2*m),inf); %输出第一层的非支配解
%     realIndividual = Archive(find(FrontNo==1),1:D);
%     realIndividual = unique(realIndividual,'rows');
    %%  Archive m+1
   
%     [FrontNo,~] = NDSort(Archive(:,D+m+1:D+2*m),inf); %输出第一层的非支配解
%     EI = Archive(find(FrontNo==1),D+m+1:D+2*m);
%      PopulationFirst = Archive(find(FrontNo==1),1:D);
%      EInum = size(EI,1)
%      EI_mid = zeros(1,m);
%      if EInum > m+1
%          for k = 1:m
%              [~,ind] = min(EI(:,k));
%              EI_mid(1,k) = (1/2)*max(EI(:,k));
% %              EI_mid(1,k) = (1/2)*(max(EI(:,k))-min(EI(:,k)));
%              realIndividual(k,:) = PopulationFirst(ind,:);
%          end
%          [~,index] = min(dist(EI,EI_mid'));
%          realIndividual(m+1,:) = PopulationFirst(index,:);
%      else
%          realIndividual = PopulationFirst;
%      end
%      realIndividual = unique(realIndividual,'rows');
   %%  first front(fitness or EI objective space)
%     [FrontNo,MaxFNo] = NDSort(EIP,inf); %输出第一层的非支配解
%     realIndividual = Population(find(FrontNo==1),:);
%     realIndividual = unique(realIndividual,'rows');
   %% fitness+EI first front
%     [FrontNum,~] = NDSort(Population,inf); %输出第一层的非支配解
%     Populationfirst = Population(find(FrontNum==1),:);
%     EIPfirst = EIP(find(FrontNum==1),:);
%     [FrontNo,~] = NDSort(EIPfirst,inf); %输出第一层的非支配解
%     realIndividual = Populationfirst(find(FrontNo==1),:);
%     realIndividual = unique(realIndividual,'rows');
     %% m+1 
     % fitness+EI first front
     [FrontNum,~] = NDSort(Population,inf); %输出第一层的非支配解
     Populationfirst = Population(find(FrontNum==1),:);
     EIPfirst = EIP(find(FrontNum==1),:);
     [FrontNo,~] = NDSort(EIPfirst,inf); %输出第一层的非支配解
     PopulationFirst = Populationfirst(find(FrontNo==1),:);
     EI = EIPfirst(find(FrontNo==1),:);

     % EI first front (fitness or EI objective space)
%      [FrontNo,~] = NDSort(EIP,inf); %输出第一层的非支配解
%      PopulationFirst = Population(find(FrontNo==1),:);
%      EI = EIP(find(FrontNo==1),:);
     
     EInum = size(EI,1)
     EI_mid = zeros(1,m);
     if EInum > m+1
         for k = 1:m
             [~,ind] = min(EI(:,k));
             EI_mid(1,k) = (1/2)*max(EI(:,k));
%              EI_mid(1,k) = (1/2)*(max(EI(:,k))-min(EI(:,k)));
             realIndividual(k,:) = PopulationFirst(ind,:);
         end
         [~,index] = min(dist(EI,EI_mid'));
         realIndividual(m+1,:) = PopulationFirst(index,:);
     else
         realIndividual = PopulationFirst;
     end
     realIndividual = unique(realIndividual,'rows');

end