function k = kernel(x, y, i, var)
    if i == 1
        k = exp((-norm(x-y)^2)/(2*var^2));
    end
    if i == 2
        k = (sum(x.*y)+1)^var
    end
end
function [eigenvalue, eigenvectors, project_invectors] = kpca(x, sigma, cls, target_dim)
    % kpca进行数据提取的函数
    psize=size(x);
    m=psize(1);     % 样本数
    n=psize(2);     % 样本维数


    % 计算核矩阵k
    l=ones(m,m);
    for i=1:m
        for j=1:m
           k(i,j)=kernel(x(i,:),x(j,:),cls,sigma); 
        end
    end


    % 计算中心化后的核矩阵
    kl=k-l*k/m-k*l/m+l*k*l/(m*m);  


    % 计算特征值与特征向量
    [v,e] = eig(kl); 
    e = diag(e);


    % 筛选特征值与特征向量
    [dump, index] = sort(e, 'descend');
    e = e(index);
    v = v(:, index);
    rank = 0;
    for i = 1 : size(v, 2)
        if e(i) < 1e-6
            break;
        else
            v(:, i) = v(:, i) ./ sqrt(e(i));
        end
        rank = rank + 1;
    end
    eigenvectors = v(:, 1 : target_dim);
    eigenvalue = e(1 : target_dim);


    % 投影
    project_invectors = kl*eigenvectors;   %计算在特征空间向量上的投影 
end

function compare
    clear all;
    close all;
    clc;
    
    % 生成非线性可分的三类数据
    if exist('X1.mat')
        load 'X1.mat'
        load 'X2.mat'
        load 'X3.mat'
        figure(1)       
        plot(X1(1, :),X1(2, :) ,'ro')
        hold on;
        plot(X2(1, :),X2(2, :),'g*')
        hold on;
        plot(X3(1, :),X3(2, :),'b.')
        hold on;
        
        title('原始数据');
        xlabel('第一维');
        ylabel('第二维');
        saveas(gcf, '原始数据图.jpg')
    else
        [X1, X2, X3] = generate_data();
        save 'X1.mat'  X1
        save 'X2.mat'  X2
        save 'X3.mat'  X3
    end

    X = [X1 X2 X3];
    [nFea, nSmps] = size(X);
    nClsSmps = nSmps / 3;
    
    % PCA
    [vec_pca, Y_pca, value_pca] = pca(X');
    Y_pca = Y_pca';
    
    figure(2);
    plot(Y_pca(1, 1 : nClsSmps),Y_pca(2, 1 : nClsSmps), 'ro');
    hold on;
    plot(Y_pca(1, nClsSmps + 1 : 2 * nClsSmps),Y_pca(2, nClsSmps + 1 : 2 * nClsSmps), 'g*');
    hold on;
    plot(Y_pca(1, 2 * nClsSmps + 1 : end),Y_pca(2, 2 * nClsSmps + 1 : end), 'b.');
    hold on;
    title('PCA');
    xlabel('第一维');
    ylabel('第二维');
    saveas(gcf, 'PCA投影图.jpg')
    
    % KPCA
    percent = 1;
    var   = 2; % 1 代表高斯核，2代表多项式核，3代表线性核
    sigma = 6; % 核参数
    [vec_KPCA, value_KPCA, Y_pca] = kpca(X', sigma, var, 2);
    Y_pca = Y_pca';
    
    figure(3);
    plot(Y_pca(1, 1 : nClsSmps),Y_pca(2, 1 : nClsSmps), 'ro');
    hold on;
    plot(Y_pca(1, nClsSmps + 1 : 2 * nClsSmps),Y_pca(2, nClsSmps + 1 : 2 * nClsSmps), 'g*');
    hold on;
    plot(Y_pca(1, 2 * nClsSmps + 1 : end),Y_pca(2, 2 * nClsSmps + 1 : end), 'b.');
    hold on;
    str = strcat('KPCA', '(p =', num2str(sigma), ')');
    title(str);
    xlabel('第一维');
    ylabel('第二维');
    str = strcat(str, '.jpg')
    saveas(gcf, str)
end