function Offspring = P_generator(MatingPool,Boundary,Coding,MaxOffspring)
% This function includes the SBX crossover operator and the polynomial
% mutatoion operator. 


    [N,D] = size(MatingPool);
    if nargin < 4 || MaxOffspring < 1 || MaxOffspring > N
        MaxOffspring = N;
    end
    
    switch Coding
        case 'Real'
            ProC = 1;       
            ProM = 1;     
            DisC = 20;   	
            DisM = 20;   	
            Offspring = zeros(N,D);
           %% SBX
            %产生N个子代
            for i = 1 : 2 : N
                beta = zeros(1,D);
                miu  = rand(1,D);
                beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
                beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
                beta = beta.*(-1).^randi([0,1],1,D);
                beta(rand(1,D)>ProC) = 1;
                Offspring(i,:)   = (MatingPool(i,:)+MatingPool(i+1,:))/2+beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
                Offspring(i+1,:) = (MatingPool(i,:)+MatingPool(i+1,:))/2-beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
            end
            %取出前N-1个子代个体作为当前代子代
            Offspring = Offspring(1:MaxOffspring,:);
            
            %% 边界设置
            if MaxOffspring == 1
                MaxValue = Boundary(1,:);%Boundary是2*D的矩阵，第一行全是1，第二行全是0
                MinValue = Boundary(2,:);
            else
                MaxValue = repmat(Boundary(1,:),MaxOffspring,1);%矩阵重复，以Boundary(1,:)为矩阵元素，重复MaxOffspring行1列，即最后的矩阵为MaxOffspring*D维的全1矩阵
                MinValue = repmat(Boundary(2,:),MaxOffspring,1);
            end
            %%
            k    = rand(MaxOffspring,D);
            miu  = rand(MaxOffspring,D);
            %多项式变异（子代所有粒子的每个维度进行变异）
            Temp = k<=ProM & miu<0.5;
            Offspring(Temp) =  Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
            Temp = k<=ProM & miu>=0.5; 
            Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
          
            %边界约束
            Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
            Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);

    end
end