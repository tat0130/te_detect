clc
clear 
close all
%% 产生训练数据
Train = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d00.dat');
Train = Train';
xtrain=Train(1:480,[1:22,42:52]);
%% 产生测试数据
Test = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d07_te.dat');
xtest =Test(1:960,[1:22,42:52]);
%标准化处理：
x_mean = mean(xtrain);                             
x_std = std(xtrain);                        
[x_row,x_col] = size(xtrain);                                                       
xtrain=(xtrain-repmat(x_mean,x_row,1))./repmat(x_std,x_row,1);
xtest_R = size(xtest,1);
xtest=(xtest-repmat(x_mean,xtest_R,1))./repmat(x_std,xtest_R,1);
%% 
CM = (xtrain'*xtrain)/(x_row-1);
[T,lamda] = eig(CM);                                                                       
E = flipud(diag(lamda));                            
num_P = 1;                                         
while sum(E(1:num_P))/sum(E) < 0.75  
num_P = num_P +1;
end                                                 
P = T(:,x_col-num_P+1:x_col);
TT=xtrain*P;
%% 
%求T2和Q统计量
[P_row,P_row1] = size(P*P');
I = eye(P_row,P_row1);
for i = 1:xtest_R
    T2(i)=xtest(i,:)*P*pinv(lamda(x_col-num_P+1:x_col,x_col-num_P+1:x_col))*P'*xtest(i,:)';
    %T1(i)=hampel(T2(i),20,1);
    SPE(i) = xtest(i,:)*(I - P*P')*(I - P*P')'*xtest(i,:)';
    %SPE1(i) = hampel(SPE(i),20,1);
end

%for i = 1:xtest_R
%    T1(i)=hampel(T2(i),20,1);
%    SPE1(i) = hampel(SPE(i),20,1);
%end

%T1=hampel(T2,4,1);
%SPE1 = hampel(SPE,4,1);

JT=num_P*(x_row-1)*(x_row+1)*finv(0.99,num_P,x_row - num_P)/(x_row*(x_row - num_P));
for i = 1:3
    theta(i) = sum((E(num_P+1:x_col)).^i);
end
h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2);
ca = norminv(0.95,0,1);
JQ = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0);                           



%T2故障检测准确率
  k1=0;
for i=160:960
   if T2(i)>JT
       k1=k1+1;
   end
end
FDRT=k1/(960-160)
%T2误报率
k1=0;
for i=1:159
   if T2(i)>JT
       k1=k1+1;
   end
end
FART=k1/159
%SPE故障检测准确率
k1=0;
for i=160:960
   if SPE(i)>JQ
       k1=k1+1;
   end
end
FDRSPE=k1/(960-160)
%SPE误报率
k1=0;
for i=1:159
   if SPE(i)>JQ
       k1=k1+1;
   end
end
FARSPE=k1/159

%% 
%给出仿真图
figure('color',[1 1 1]);
subplot(2,1,1)   
plot(T2,'linewidth',1.5,'color','b');
hold on
plot(repmat(JT,1,xtest_R),'r-.','linewidth',1); 
xlabel('Sample number');
ylabel('T^2');
legend('Statistic','Control limit');
title(['检测准确率为',num2str(FDRT)],['误报率为',num2str(FART)]);

subplot(2,1,2)   
plot(SPE,'linewidth',1.5,'color','b');
hold on
plot(repmat(JQ,1,xtest_R),'r-.','linewidth',1); 
xlabel('Sample number');
ylabel('SPE');
legend('Statistic','Control limit');
title(['检测准确率为',num2str(FDRSPE)],['误报率为',num2str(FARSPE)]);

