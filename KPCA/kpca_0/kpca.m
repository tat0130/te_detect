
    clear all;
    close all;
    clc;
    
    % ���ɷ����Կɷֵ���������
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
        
        title('ԭʼ����');
        xlabel('��һά');
        ylabel('�ڶ�ά');
        saveas(gcf, 'ԭʼ����ͼ.jpg')
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
    xlabel('��һά');
    ylabel('�ڶ�ά');
    saveas(gcf, 'PCAͶӰͼ.jpg')
    
    % KPCA
    percent = 1;
    var   = 2; % 1 �����˹�ˣ�2�������ʽ�ˣ�3�������Ժ�
    sigma = 6; % �˲���
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
    xlabel('��һά');
    ylabel('�ڶ�ά');
    str = strcat(str, '.jpg')
    saveas(gcf, str)
