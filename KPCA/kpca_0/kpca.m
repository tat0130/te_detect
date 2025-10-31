
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
