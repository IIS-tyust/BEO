function  [Population,prefitnessP,realIndividual] = KDPEA_LCBEI(dmodel,Population,prefitnessP,Boundary,reference_vector,Archive_train)
% 返回下一代父本和真实评价的个体
%%
%-------------------参数设置----------------
format compact;
Coding = 'Real';
[N,m] = size(reference_vector);
D = size(Population,2);
ER = [];

%------------------找每个参考向量上对应的参考点(理想点)-------
ref_point = reference_point_assign(prefitnessP,reference_vector);
%------------------SBX and PM 产生候选解---------------
 [MatingPool] = F_mating(Population);
%     MatingPool= F_mating_new(Population,N);
Offspring = P_generator(MatingPool,Boundary,Coding,N);
prefitnessO = zeros(size(Offspring,1),m);
ER = prefitnessO; %ER均方误差
%---------使用Kriling模型估计子代(候选解)-----------
for i = 1:size(Offspring,1)
    for k = 1:m
        [prefitnessO(i,k),~,~] = predictor(Offspring(i,:),dmodel(k));
    end
end
Of = [Offspring prefitnessO];
Xf = [Population prefitnessP];
f = [prefitnessP;prefitnessO];
XO = [Xf;Of];
%----------------按角度关联参考向量------------
indiv_ind = Angle_assign_refvec(f,reference_vector);%每个参考向量可能有多个个体关联
%------------------产生新参考点-----------------
ref_point_new = reference_point_assign(prefitnessO,reference_vector);
%--------------------更新每个参考向量上的参考点------
ref_point = selection_reference_point(ref_point,ref_point_new);
%         % -------------------更新每个参考向量上的理想点--------------
%         ideal_point = ideal_point_generation(indiv_ind,f,reference_vector);
%--------------------按照指标选择个体-------------
%     Xf = selection_individual(indiv_ind,XO,reference_vector,ideal_point);
Xf = [];
Xf = selection_individual(indiv_ind,XO,reference_vector,ref_point);
Population = [];
prefitnessP = [];
Population = Xf(:,1:D); %根据指标选择出来的个体作为下一代父本
prefitnessP = Xf(:,D+1:D+m);
%--------组建EI矩阵，根据EI值排序，选择较大的1/2个体进行真实评价
f_min = min(Archive_train(:,D+1:D+m),[],1);
% EI_times = [];
% prefitness = [];
for i = 1:size(Population,1)
    for t = 1:m
        [prefitness(i,t),~,ER(i,t)] = predictor(Population(i,:),dmodel(t));
        s(i,t) = sqrt(ER(i,t));
        EI(i,t) = (f_min(1,t) - prefitness(i,t)) * Gaussian_CDF((f_min(1,t)-prefitness(i,t))./ s(i,t))...
            + s(i,t) * Gaussian_PDF((f_min(1,t)-prefitness(i,t))./ s(i,t));
        LCB(i,t) = LCB_function(prefitness(i,t),ER(i,t));
    end
    EI_times(i,1) = prod(EI(i,:));
    LCB_times(i,1) = prod(LCB(i,:));
end
% [EI_value,EI_id] = sort(EI_times,'descend');%EI 选取较大
% [minlcb,minlcbid] = sort(LCB_times, 'ascend');%LCB选取较小
EILCB = [Population -EI_times LCB_times];
[nospop,rank_info] = non_dominant_sort(EILCB,2,D);
realIndividual = nospop(1:rank_info(1,1),1:D);
% num = round((1/2)* size(EI_id,2));
% num = 5;
% realIndividual = Population(EI_id(1,1:num),:);
end





