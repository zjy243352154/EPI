q = 1;
    clear pyof_mat;
    par_vec = par_mat(q,:);

L1 = zeros(4*870*2,870*5);%四年的部门产出约束 4*870*2
L2 = zeros(5*29*2,870*5);%五年的行业产出约束 5*29*2
L3 = zeros(5*30*1,870*5);%五年的省份产出约束 5*30*1
L4 = zeros(5*1,870*5);%五年的GDP增长约束   5*1
L5 = zeros(5*870,870*5);%五年的投入产出平衡约束 5*870
L6 = zeros(5*31*1,870*5); %五年的能耗约束 5*31*1
L7 = zeros(5*1,870*5);%NOx排放量总量约束 5*1
L8 = zeros(5*1,870*5);%COD排放量总量约束 5*1
L9 = zeros(5*1,870*5);%AN排放量总量约束 5*1
L10 = zeros(5*1,870*5);%用水量总量约束 5*1
L11 = zeros(5*1,870*5);% SO2总量约束 5*1
L12 = zeros(5*1,870*5);% SD总量约束 5*1
L13 = zeros(1*870*2,870*5);
% RHS
R1=zeros(4*870*2,1);%四年的部门约束 5*870*2
R2=zeros(5*29*2,1);%五年的行业产出约束 5*29*2
R3=zeros(5*30*1,1);%五年的省份产出约束 5*30*2
R4=zeros(5*1,1);%五年的GDP增长率约束   5*1
R5=zeros(5*870,1);%五年的投入产出平衡约束5*870
R6 = zeros(5*31,1);%五年的能耗约束5*30
R7 = zeros(5*1,1);
R8 = zeros(5*1,1);
R9 = zeros(5*1,1);
R10 = zeros(5*1,1);
R11 =  zeros(5*1,1);
R12 =  zeros(5*1,1);
R13 = zeros(870*2,1);
for i = 1:5
    matrix_temp{i} = reshape([870*i-869:870*i],[29,30]);
end

%%写入约束矩阵
%%部门产出约束
%%L1
%产出下限
for i =1:4
    for l = 1:870
        [s,r] = find(matrix_temp{1} == l);
        if s == 1 || s==28 || s == 29
            L1((i-1)*870+l,i*870+l)=-1;
            L1((i-1)*870+l,(i-1)*870+l) = par_vec(3);
        else
            L1((i-1)*870+l,i*870+l)=-1;
            L1((i-1)*870+l,(i-1)*870+l) = par_vec(1);
        end
    end
end


%产出上限
for i= 1:4
    for l = 1:870
        [s,r] = find(matrix_temp{1}== l);
        if s == 1 || s==28 || s == 29
            L1((3+i)*870+l,i*870+l) =1;
            L1((3+i)*870+l,(i-1)*870+l) = -par_vec(4);%产出下限
        else
            L1((3+i)*870+l,i*870+l) =1;
            L1((3+i)*870+l,(i-1)*870+l) = -par_vec(2);%产出下限
        end
    end
end

%%R1
R1=R1;


%%行业产出约束
%L2
%行业产出下限

for r = 1:29
    L2(r,matrix_temp{1}(r,:))=-1;
end

for i = 1:4
    for r = 1:29
        L2(i*29+r,matrix_temp{i+1}(r,:)) = -1;
        L2(i*29+r,matrix_temp{i}(r,:)) = par_vec(5);
    end
end

for r = 29*5+1:29*6
    L2(r,matrix_temp{1}(r-29*5,:))=1;
end

for i = 1:4
    for r = 1:29
        L2((5+i)*29+r,matrix_temp{i+1}(r,:))=1;
        L2((5+i)*29+r,matrix_temp{i}(r,:))=-par_vec(6);
    end
end

%R2 0.864/1.168
for r = 1:29
    R2(r) = sum(X_2025(matrix_temp{1}(r,:)))*(-par_vec(5));
    R2(r+29*5) = sum(X_2025(matrix_temp{1}(r,:)))*par_vec(6);
end

%%省份产出约束

for s = 1:30
    L3(s,matrix_temp{1}(:,s))=-1;
end
for i = 1:4
    for s= 1:30
        L3(i*30+s,matrix_temp{i+1}(:,s)) = -1;
        L3(i*30+s,matrix_temp{i}(:,s)) = par_vec(7);
    end
end



for s = 1:30
    R3(s) = sum(X_2025(matrix_temp{1}(:,s)))*(-par_vec(7));
end

%%GDP增长率约束

L4(1,matrix_temp{1}(:)) = -B_index;
for i = 1:4
    L4(i+1,matrix_temp{i+1}(:)) = -B_index;
    L4(i+1,matrix_temp{i}(:)) = par_vec(8)*B_index;
end
R4(1) = (B_index)*(-par_vec(8))*X_2025;

%%投入产出平衡约束
%
% Amat2 = Amat-diag(ones(1,870));

% R5 = -repmat(F_2025,5,1);
Amat3 = Amat - diag(ones(1,870));
temp3 = find((-Amat3)*X_2025 < F_2025 );
Amat3(temp3,:) = 0;
L5 = [];
for i = 1:5
    L5 = blkdiag(L5,Amat3);
end
F_2025(temp3) = 0;
R5 = -repmat(F_2025,5,1);


%%能耗约束
energy_uplmt = [568000:8000:600000];% 设置全国能耗上限
for i = 1:5
    for j = 1:30
        L6((i-1)*30+j,matrix_temp{i}(:,j)) = EPI{i}(8,matrix_temp{1}(:,j));
        R6((i-1)*30+j) = energy_uplmt(i)*scale(j)*0.9*1.5;
    end
end
for i = 1:5
    L6(150+i,matrix_temp{i}) = EPI{i}(8,:);
    R6(150+i) = energy_uplmt(i)*0.9;
end

% NOX约束
nox_uplmt = linspace(paifang2025(3),0.9*paifang2025(3),5);
for i = 1:5
    L7(i,matrix_temp{i}(:)) = EPI{i}(3,matrix_temp{1}(:));
    R7(i) = nox_uplmt(i);
end
% COD约束
COD_uplmt = linspace(paifang2025(5),paifang2025(5)*0.92,5);
for i = 1:5
    L8(i,matrix_temp{i}(:)) = EPI{i}(5,matrix_temp{1}(:));
    R8(i) = COD_uplmt(i);
end
% AN约束
an_uplmt = linspace(paifang2025(6),paifang2025(6)*0.92,5);
for i = 1:5
    L9(i,matrix_temp{i}(:)) = EPI{i}(6,matrix_temp{1}(:));
    R9(i) = an_uplmt(i);
end
% SO2约束
so2_uplmt = linspace(paifang2025(2),paifang2025(2),5);
for i = 1:5
    L10(i,matrix_temp{i}(:)) = EPI{i}(2,matrix_temp{1}(:));
    R10(i) = so2_uplmt(i);
end
% Sd约束
sd_uplmt = linspace(paifang2025(4),paifang2025(4),5);
for i = 1:5
    L11(i,matrix_temp{i}(:)) = EPI{i}(4,matrix_temp{1}(:));
    R11(i) = sd_uplmt(i);
end
%%水资源约束
wat_uplmt = linspace(paifang2025(1),paifang2025(1),5);
for i = 1:5
    L12(i,matrix_temp{i}(:)) = EPI{i}(1,matrix_temp{1}(:));
    R12(i) = wat_uplmt(i);
    
end

for l =1:870
    [s,r] = find(matrix_temp{1} == l);
    if s == 1 || s==28 || s == 29
        ub_temp(l) = X_2025(l)*1.1;
        lb_temp(l) = X_2025(l)*0.9;
    else
        ub_temp(l) = X_2025(l)*1.2;
        lb_temp(l) = X_2025(l)*0.8;
    end
end
ub = [ub_temp';inf(870*4,1)];
lb = [lb_temp';zeros(870*4,1)];
%%
for i = 1:870
L13(i,i) = 1;
L13(870+i,i) = -1;
R13(i) = ub(i);
R13(i+870) = lb(i);
end

%%总约束


LHS = [L1;L2;L3;L4;L5;L6;L7;L8;L9;L10;L11;L12];
RHS = [R1;R2;R3;R4;R5;R6;R7;R8;R9;R10;R11;R12];

%%
LHS_temp = [L1;L2;L3;L4;L6;L7;L8;L9;L10;L11;L12];
RHS_temp = [R1;R2;R3;R4;R6;R7;R8;R9;R10;R11;R12];
%% 

%%目标系数
% for j = 1:8
%     C{j} = zeros(1,5*870);
%     for i = 1:5
%         C{j}(1,870*i-869:870*i) = EPI{i}(j,:);
%     end
% end
% for i = 1:5
%     C{9}(1,870*i-869:870*i) = -B_index;
% end
%% 求解单目标
options = optimoptions('linprog','Algorithm','interior-point','MaxIterations',Inf,'OptimalityTolerance',1.0000e-10);
tic
for j = 1:9
    if j <= 8
        c = C{j};
        LHS_limit = LHS;
        RHS_limit = RHS;
        limit_sector = find(paifang_2017(j,:)==0);
        for i = 1:4
            for l = 1:length(limit_sector)
                LHS_limit(870*(i-1)+limit_sector(l),(i-1)*870+limit_sector(l)) = 0.999;
                LHS_limit(870*(i-1)+limit_sector(l),i*870+limit_sector(l)) = -1;
                LHS_limit(870*(i+3)+limit_sector(l),i*870+limit_sector(l)) =  1;
                LHS_limit(870*(i+3)+limit_sector(l),(i-1)*870+limit_sector(l)) =  -1.001;
            end
        end
        A = LHS_limit;
        B = RHS_limit;
        lb_limit =lb;
        ub_limit = ub;
        lb_limit(limit_sector) =X_2025(limit_sector)*0.999;
        ub_limit(limit_sector) =X_2025(limit_sector)*1.001;
        [x,fval,exitflag,output] = linprog(c,A,B,[],[],lb_limit,ub_limit,options);
        output
        result1{j}{1} = x;
        j
    else
        c = C{j};
        A = LHS;
        B =RHS;
        [x,fval,exitflag,output] = linprog(c,A,B,[],[],lb',ub',options);
        exitflagg{j} = exitflag;
        output
        result1{j}{1} = x;
        j
    end
end
toc
for i = 1:9
    xlswrite(['solution',num2str(q),'.xlsx'],result1{i}{1},lgname{i});
end


%% 根据LX计算支付矩阵
global pyof_mat;
for j = 1:9
    for l  = 1:5*870
        [ll,i] = find(reshape([1:6*870],870,6)==l);
        result{i}{j}{1}(ll,1) = result1{j}{1}(l);
        result{i}{10}{1} = zeros(870,1);
    end
end
EPI;%{5}(8*870)
%X{9}(870*5,1)
for j = 1:9
    for i = 1:5
        X_temp{i}{j} = result{i}{j}{1};
    end
end
%payoff_matrix
for j = 1:8
    for jj = 1:8
        EP{j}(jj)=0;
        for i = 1:5
            EP{j}(jj) =  sum(EPI{i}(jj,:)*(X_temp{i}{j}))+EP{j}(jj);
        end
        pyof_mat(j,jj) = EP{j}(jj);
    end
end
for j = 1:9
    pyof_mat(j,9) = -C{9}*result1{j}{1};
end
for j = 1:8
    pyof_mat(9,j) =0;
    pyof_mat(9,j) = C{j}*result1{9}{1};
end
xlswrite("pyof_mat.xlsx",pyof_mat,'1','A1:I9');
%% 求解多目标
start_x = result1{1}{1}';
for i = 1:8
    start_x = [start_x;result1{i+1}{1}'];
end
for i = 1:50
tic
    options19 = optimoptions("ga","Display","iter","ConstraintTolerance",0.01,"PopulationSize",9,'InitialPopulationMatrix',start_x,'FunctionTolerance',1e-5,...
        'MaxStallGenerations',10);
    % 求解
    [cpx,cpfval,exitflag,output] = ga(@ideal_point,l,LHS_temp,RHS_temp,[],[],lb,ub,[],options19);
    toc
     i
    exitflag_cell{i} = exitflag;
    comp_cell{i} = cpx;
    compfval_cell{i} = cpfval;
end