clc 
clear all
%导入数据
load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d00.dat')
Test=load('C:\Users\admin\Desktop\data set\tennessee-eastman-profBraatz-master\d21_te.dat')
% % d00=xlsread('C:\Users\PinePineSong\Desktop\future\PROGRAM\TE\训练集1');
% % d01_te=xlsread('C:\Users\PinePineSong\Desktop\future\PROGRAM\TE\测试集1');
d00_te=Test';

%主元确定
gx=0;%实际贡献率
GXset=0.85;%期望贡献率
L=1;%初始主元个数
a=0.15;%置信水平1-a
%%数据归一化
%zbarobs离线
Zbarobs=mean(d00,2);
for m=1:1:size(d00,1)
   for N=1:1:size(d00,2)
       deltaobs(m,N)=(d00(m,N)-Zbarobs(m,1));
   end 
end
%delta2obs
for ms=1:1:size(d00,1)
       delta2obs(ms,1)=deltaobs(ms,:)*deltaobs(ms,:)'/(size(d00,2)-1); 
end
for N=1:1:size(d00,2)
    Z(:,N)=(d00(:,N)-Zbarobs(:,1))./sqrt(delta2obs(:,1));
end
% % % SVD
% % [U,Gama,V]=svd(Z*Z'/(size(d00,1)-1));
%eig
[P,Gama]=eig(Z*Z'/(size(d00,1)-1));
lambda=wrev(diag(Gama));
Gama=diag(lambda);
P=fliplr(P);
for m=1:1:size(d00,1)
    delta2(m,1)=Gama(m,m);
end
%确定主元个数
while gx<=GXset
   gx=sum(lambda(1:L,:))/sum(lambda(:,:));
   L=L+1;
end
for j=L+1:1:size(d00,1)
    th1(:,j)=(delta2(j,1))^1;
    th2(:,j)=(delta2(j,1))^2;
    th3(:,j)=(delta2(j,1))^3;
end
theta1=sum(th1);
theta2=sum(th2);
theta3=sum(th3);
h0=1-2*theta1*theta3/(3*theta2*theta2);
%阈值生成
ca=norminv(1-a);
Fa=finv(1-a,L,size(d00,2)-L);
JthSPE=theta1*((ca*sqrt(2*theta2*h0*h0)/theta1)+1+(theta2*h0*(h0-1))/(theta1*theta1))^(1/h0)*ones(1,size(d00_te,2));
JthT2PCA=(L*(size(d00,2)^2-1)*Fa)/(size(d00,2)*(size(d00,2)-L))*ones(1,size(d00_te,2));
%在线测试
for i=1:1:size(d00_te,2)
    znorm(:,i)=(d00_te(:,i)-Zbarobs(:,1))./sqrt(delta2obs(:,1));
    SPE(:,i)=znorm(:,i)'*(eye(size(d00_te,1))-P(:,1:L)*P(:,1:L)')*znorm(:,i);
    T2PCA(:,i)=znorm(:,i)'*P(:,1:L)*inv(diag(lambda(1:L,:)))*P(:,1:L)'*znorm(:,i);
end

%% 图
k=1:1:size(d00_te,2);
figure(1)
subplot(2,1,1)   
plot(k,SPE(:,k),'b',k,JthSPE(:,k),'-r','LineWidth',0.2);
% axis([0 size(d00_te,2) 0 10]); %坐标范围
% set(gca,'xtick',0:100:size(d00_te,2)); %显示间隔x
% set(gca,'ytick',0:1:10); %显示间隔y
legend('SPE','JthSPE','Location','northwest' );
title('TEST--SPE');
xlabel('$k$','interpreter','latex');
ylabel('JthSPE-SPE($k$)','interpreter','latex');
hold on;   
subplot(2,1,2)   
plot(k,T2PCA(:,k),'b',k,JthT2PCA(:,k),'-r','LineWidth',0.2);
% axis([0 size(d00_te,2) 0 10]); %坐标范围
% set(gca,'xtick',0:100:size(d00_te,2)); %显示间隔x
% set(gca,'ytick',0:1:10); %显示间隔y
legend('T2PCA','JthT2PCA','Location','northwest' );
title('TEST--T2PCA');
xlabel('$k$','interpreter','latex');
ylabel('JthT2PCA-T2PCA($k$)','interpreter','latex');
hold on;  