clc
clear 
close all
%% 产生训练数据
Train = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d00.dat');
Train = Train';
xtrain=Train(1:480,[1:22,42:52]);
%% 产生测试数据
Test = load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d01_te.dat');
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
while sum(E(1:num_P))/sum(E) < 0.85  
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
    SPE(i) = xtest(i,:)*(I - P*P')*(I - P*P')'*xtest(i,:)';                                                                                    
end
                      
JT=num_P*(x_row-1)*(x_row+1)*finv(0.999,num_P,x_row - num_P)/(x_row*(x_row - num_P));
for i = 1:3
    theta(i) = sum((E(num_P+1:x_col)).^i);
end
h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2);
ca = norminv(0.999,0,1);
JQ = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0);                           

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

subplot(2,1,2)   
plot(SPE,'linewidth',1.5,'color','b');
hold on
plot(repmat(JQ,1,xtest_R),'r-.','linewidth',1); 
xlabel('Sample number');
ylabel('SPE');
legend('Statistic','Control limit');
