
function  [Population,fitness,real_individual] = K_RPEA(dmodel,Population,fitness,Boundary,reference_vector,Archive_train)
% 返回下一代父本和真实评价的个体
%%
%-------------------参数设置----------------
format compact;
Coding = 'Real';
maxFEK = 100000;
N = size(reference_vector,1);
m = size(fitness,2);
D = size(Population,2);
EP = [];%存储非支配解,维数为m
FEK = 0;
ER = zeros(size(Population,1),m);
 %------------------找每个参考向量上对应的参考点(理想点)-------
ref_point = reference_point_assign(fitness,reference_vector);
%%     % ---------------------主循环-----------------------
while FEK <= maxFEK
    %------------------SBX and PM---------------
    [MatingPool] = F_mating(Population);
%     MatingPool= F_mating_new(Population,N);
    Offspring = P_generator(MatingPool,Boundary,Coding,N);
    prefitness = zeros(size(Offspring,1),m);
    ER = prefitness; %ER均方误差
    %---------使用Kriling模型估计子代-----------
    for i = 1:size(Offspring,1)
        for k = 1:m
            [prefitness(i,k),~,ER(i,k)] = predictor(Offspring(i,:),dmodel(k));
        end
    end
    FEK = FEK + size(Offspring,1);
  
    Of = [Offspring prefitness];
    Xf = [Population fitness];
    f = [fitness;prefitness];
    XO = [Xf;Of];
    %----------------按角度关联参考向量------------
    indiv_ind = Angle_assign_refvec(f,reference_vector);%每个参考向量可能有多个个体关联
    %------------------产生新参考点-----------------
    ref_point_new = reference_point_assign(prefitness,reference_vector);
    %--------------------更新每个参考向量上的参考点------
    ref_point = selection_reference_point(ref_point,ref_point_new);
    %         % -------------------更新每个参考向量上的理想点--------------
    %         ideal_point = ideal_point_generation(indiv_ind,f,reference_vector);
    %--------------------按照指标选择个体-------------
    %     Xf = selection_individual(indiv_ind,XO,reference_vector,ideal_point);
    Xf = selection_individual(indiv_ind,XO,reference_vector,ref_point);
    Population = zeros(size(Population,1),D);
    fitness = zeros(size(Population,1),m);
    Population = Xf(:,1:D);
    fitness = Xf(:,D+1:D+m);
end
 %--------组建EI矩阵，根据EI值排序，选择较大的1/2个体进行真实评价
   f_min = min(Archive_train(:,D+1:D+m),[],1);
   for i = 1:size(Population,1)
        for k = 1:m
            [prefitness(i,k),~,ER(i,k)] = predictor(Population(i,:),dmodel(k));
            s(i,k) = sqrt(ER(i,k));
            EI(i,k) = (f_min(1,k) - prefitness(i,k)) * Gaussian_CDF((f_min(1,k)-prefitness(i,k))./ s(i,k))...
                + s(i,k) * Gaussian_PDF((f_min(1,k)-prefitness(i,k))./ s(i,k));
        end
        EI_times(i) = prod(EI(i,:));
    end
   [EI_value,EI_id] = sort(EI_times,'descend');
    num = round((1/2)* size(EI_id,2));
    real_individual = [];
    real_individual = Population(EI_id(1,1:num),:);
end





