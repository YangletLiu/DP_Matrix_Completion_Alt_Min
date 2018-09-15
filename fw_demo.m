filename = 'H:\GitHub\DP_Matrix_Completion_Alt_Min\data\input.xlsx'; %数据文件名

P(:,:)=xlsread(filename,1,'B2:OK6041');%读取数据文件

[m n] = size(P); %返回data数据文件里的矩阵大小
fprintf('data has been loaded: m = %d, n = %d; \n', m,n);


