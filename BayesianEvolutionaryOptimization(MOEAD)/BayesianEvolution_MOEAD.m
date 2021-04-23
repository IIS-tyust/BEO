%%
clear;
clc;
addpath(genpath('dace'));
clear dmodel;
warning off;

FEmax = 300;
m = 3;
N = 50; 

for tt = 1:1
    tt 
    switch tt
        case 1
            problem = 'DTLZ1';
%              problem = 'ZDT1';
        case 2
%             problem = 'DTLZ2';
             problem = 'ZDT2';
        case 3
%             problem = 'DTLZ3';
             problem = 'ZDT3';
        case 4
%             problem = 'DTLZ4';
             problem = 'ZDT4';
        case 5
%             problem = 'DTLZ5';
             problem = 'ZDT6';
%         case 6
%             problem = 'DTLZ6';
%         case 7
%             problem = 'DTLZ7';
    end
    for mmm =1:1
        mmm
      %% -----------初始化阶段--------
%         [reference_vector,N] = UniformPoint(N,m);
%         reference_vector(reference_vector == 0) = 0.000001;
        % N = size(reference_vector(:,1),1);
        [~,Boundary,Coding] = P_objective('init',problem,m,N);
        D = size(Boundary,2);
        sample_size = 11*D-1; %训练样本大小
        Xn = lhsdesign(sample_size,D); %拉丁超立方体初始化
        ub = Boundary(1,:);%决策变量所有维的上界
        lb = Boundary(2,:);
        Population = bsxfun(@plus,lb,bsxfun(@times,Xn,(ub-lb)));%产生训练样本
        fitness = P_objective('value',problem,m,Population);
        FE = size(fitness,1);
        Archive_train = [Population fitness];
        Archive_train_Dec = Archive_train(:,1:D);
%         Archive_train = unique(Archive_train,'rows');
       [~,distinct] = unique(Archive_train_Dec,'rows');%distinct返回population中的个体(按每维升序排序)所在种群序号
       Archive_train_new = Archive_train(distinct,:);
        [FrontNo,MaxFNo] = NDSort(Archive_train_new(:,D+1:D+m),inf); %输出第一层的非支配解
        Archive_data_new = Archive_train_new(find(FrontNo==1),:); 
    
        %-------------训练模型------------
        theta1 = 0.5*ones(1,D); lob = 0.00001*ones(1,D); upb = 200*ones(1,D);
        % theta1 = 10*ones(1,D); lob = 0.000001*ones(1,D); upb = 20*ones(1,D);
        for t = 1:m
            dmodel(t) = dacefit(Archive_train_new(:,1:D),Archive_train_new(:,D+t), @regpoly0, @corrgauss, theta1, lob, upb);
        end
      %% ----------Evolutionary Initioanlization-----------
    
       %% Generate the weight vectors
       [W,N] = UniformPoint(N,m);
        T = ceil(N/10);

      %% Detect the neighbours of each solution
       B = pdist2(W,W); %计算W之间的欧氏距离
       [~,B] = sort(B,2); %B中的各行元素按升序排序
       B = B(:,1:T); %取每行最小的T个距离值的W索引号
    
     %% Generate random population
%       [Population,Boundary,~] = P_objective('init',problem,m,N);
       RandList = randperm(sample_size);
       Population = Population(RandList(1:N),:);
      %模型估计
      fitnessP = zeros(N,m); ER = fitnessP; s = ER; EIP = fitnessP;
      f_min = min(Archive_train(:,D+1:D+m),[],1);
      for i = 1:N
          for t = 1:m
              [fitnessP(i,t),~,ER(i,t)] = predictor(Population(i,:),dmodel(t));
              s(i,t) = sqrt(ER(i,t));
              EIP(i,t) = (f_min(1,t) - fitnessP(i,t)) * Gaussian_CDF((f_min(1,t)-fitnessP(i,t))./ s(i,t))...
                  + s(i,t) * Gaussian_PDF((f_min(1,t)-fitnessP(i,t))./ s(i,t));
          end
      end
   %%  fitness and EI space
      Z = min(fitnessP,[],1);
%       Z = min(EIP,[],1);
        %% Main loop
        while FE <= FEmax
            FE
            %--------种群进化和产生真实评价个体------
            [realIndividual,Population,fitnessP,EIP,Z] = evolve_MOEAD(dmodel,Archive_train_new,Population,fitnessP,EIP,Boundary,W,B,T,Z);
            realnum = size(realIndividual,1)
            f_real = P_objective('value',problem,m,realIndividual);
            
            %用于测试
            %      for i = 1:size(f_real,1)
            %          for k = 1:m
            %              [prefitness(i,k),~,~] = predictor(realIndividual(i,:),dmodel(k));
            %          end
            %      end
            
            FE = FE + size(f_real,1);
            %--------------更新训练集-----------
            Archive_train = [Archive_train_new; realIndividual f_real];
             Archive_train_Dec = Archive_train(:,1:D);
%            Archive_train_size1 = size(Archive_train,1)
%             Archive_train = unique(Archive_train,'rows');
            [~,distinct] = unique(Archive_train_Dec,'rows');%distinct返回population中的个体(按每维升序排序)所在种群序号
            Archive_train_new = Archive_train(distinct,:);
%            Archive_train_size = size(Archive_train_new,1)
            %---------------更新模型------------
            for t = 1:m
                dmodel(t) = dacefit(Archive_train_new(:,1:D),Archive_train_new(:,D+t), @regpoly0, @corrgauss, theta1, lob, upb);
            end
            %-------------更新非支配解集---------
            Archive_data = [Archive_data_new; realIndividual f_real];
            Archive_data_Dec = Archive_data(:,1:D);
%             Archive_data = unique(Archive_data,'rows');
            [~,distinct] = unique(Archive_data_Dec,'rows');%distinct返回population中的个体(按每维升序排序)所在种群序号
            Archive_data_new = Archive_data(distinct,:);
            [FrontNo,MaxFNo] = NDSort(Archive_data_new(:,D+1:D+m),inf); %输出第一层的非支配解
            Archive_data_new = Archive_data_new(find(FrontNo==1),:); 
            %%  画图
%                subplot(2,2,1);
%                plot3(prefitness(:,D+1),prefitness(:,D+2),prefitness(:,D+3),'bo');
%                subplot(2,2,2);
%                plot3(f_real(:,D+1),f_real(:,D+2),f_real(:,D+3),'ro');
%             
                 PF = P_objective('true',problem,m,10000);
             
%                  subplot(2,2,3);
%                  plot3(Archive_data(:,D+1),Archive_data(:,D+2),Archive_data(:,D+3),'go');
%                  xlabel('f_1');
%                  set(gca,'XDir','reverse');
%                  ylabel('f_2');
%                  set(gca,'YDir','reverse');
%                  zlabel('f_3');
%                  title(problem);
            
%                  subplot(2,2,4);
                 IGD_Score(tt,1) = IGD(Archive_data_new(:,D+1:D+m),PF)
                 plot(FE,IGD_Score,'b.');
                 hold on
                 xlabel('FE');
                 ylabel('IGD');
                 title(problem);
                 M(FE) = getframe;
        end
        %%  用于测试    
%         PF = P_objective('true', problem,m,10000);
%         IGD_Score(tt,mmm) = IGD(Archive_data_new(:,D+1:D+m),PF);
    end
%     IGD_mean(tt,1) = mean(IGD_Score(tt,:))
%     IGD_std(tt,1)  = std(IGD_Score(tt,:));
%     save('IGD_result1.mat','IGD_Score','IGD_mean','IGD_std');
end
