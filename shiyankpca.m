tic
clc 
clear
close all
%% 产生训练数据
Train = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d00.dat');
Train = Train';
Test = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d21_te.dat');
xtrain=Train(1:480,[1:22,42:52]);
xtest =Test(1:960,[1:22,42:52]);
xtest=hampel(xtest,3,1);
xtrain1=xtrain;                         
[xtrain_row,xtrain_col] = size(xtrain);                                                      
xtrain=(xtrain-repmat(mean(xtrain1),xtrain_row,1))./repmat(std(xtrain1),xtrain_row,1);
xtest_row = size(xtest,1);
xtest=(xtest-repmat(mean(xtrain1),xtest_row,1))./repmat(std(xtrain1),xtest_row,1);

%% 计算核矩阵 
sita=10000;
%sita=xtrain_row*xtrain_col;
for i=1:xtrain_row
   for j=1:xtrain_row
  KX(i,j) = exp(-(xtrain(i,:)-xtrain(j,:))*(xtrain(i,:)-xtrain(j,:))'/sita);
   end
end
%% 中心化核矩阵
 n1= ones(xtrain_row,xtrain_row)/xtrain_row;
KX1=KX;  
KX=KX-n1*KX-KX*n1+n1*KX*n1;

for i=1:xtest_row
   for j=1:xtrain_row
  KXnew(i,j) = exp(-(xtest(i,:)-xtrain(j,:))*(xtest(i,:)-xtrain(j,:))'/sita);
   end
end
%% 中心化处理
 q1= ones(xtest_row,xtrain_row)/xtrain_row;
KXnew=KXnew-q1*KX1-KXnew*n1+q1*KX1*n1;
%% 协方差矩阵的特征值分解 
CM = KX/(xtrain_row-1);
[T,lamda] = eig(CM);                                                                       
E = flipud(diag(lamda));                            
%% 确定主元个数
num_P = 1;                                         
while sum(E(1:num_P))/sum(E) < 0.95  
num_P = num_P +1;
end                                                 
P = T(:,xtrain_row-num_P+1:xtrain_row);
%% 计算T2及SPE以及控制限
temp1=KXnew*T; temp2=KXnew*T(:,xtrain_row-num_P+1:xtrain_row);
for i = 1:xtest_row
    T2(i)=KXnew(i,:)*P*pinv(lamda(xtrain_row-num_P+1:xtrain_row,xtrain_row-num_P+1:xtrain_row))*P'*KXnew(i,:)';  
    SPE(i) = temp1(i,:)*temp1(i,:)'-temp2(i,:)*temp2(i,:)';                                                                                    
end    

%T1=T2;
%SPE1=SPE;
%T1=hampel(T2,10,0.01);
%SPE1 = hampel(SPE,10,0.01);
%T1=smoothts(T2);
%SPE1 = smoothts(SPE);

temp1_obs=KX*T; temp2_obs=KX*T(:,xtrain_row-num_P+1:xtrain_row);
for i = 1:xtrain_row
    T2_obs(i)=KX(i,:)*P*pinv(lamda(xtrain_row-num_P+1:xtrain_row,xtrain_row-num_P+1:xtrain_row))*P'*KX(i,:)';  
    SPE_obs(i) = temp1_obs(i,:)*temp1_obs(i,:)'-temp2_obs(i,:)*temp2_obs(i,:)';                                                                                    
end
JT=ksdensity(T2_obs,0.95,'function','icdf');
JQ=ksdensity(SPE_obs,0.95,'function','icdf');                         


%% 
% % 计算故障检测率和故障误报率
%故障检测率
FDR1=0;
for i=160:960
    if JT < T2(i) 
        FDR1 = FDR1+1;
    end
end
FDR1=FDR1/(960-160);

%故障误报率
FAR1=0;
for i=1:160  
    if JT< T2(i) 
        FAR1 = FAR1+1;
    end
end
FAR1=FAR1/160;


%故障检测率
FDR2=0;
for i=160:960
    if JQ < SPE(i) 
        FDR2 = FDR2+1;
    end
end
FDR2=FDR2/(960-160);

%故障误报率
FAR2=0;
for i=1:160  
    if JQ< SPE(i) 
        FAR2 = FAR2+1;
    end
end
FAR2=FAR2/160;

%% 
figure('color',[1 1 1]);
subplot(2,1,1)   
plot(T2,'linewidth',1.5,'color','b');
hold on
plot(repmat(JT,1,xtest_row),'r-.','linewidth',1);  
xlabel('Sample number');
ylabel('T^2');
legend('Statistic','Control limit');
title(['检测准确率为',num2str(FDR1)],['误报率为',num2str(FAR1)]);

subplot(2,1,2)   
plot(SPE,'linewidth',1.5,'color','b');
hold on
plot(repmat(JQ,1,xtest_row),'r-.','linewidth',1);  
xlabel('Sample number');
ylabel('SPE');
legend('Statistic','Control limit');
title(['检测准确率为',num2str(FDR2)],['误报率为',num2str(FAR2)]);

toc