% 加载.mat文件（如果尚未加载）
load('Yale.mat');

% 查看X的结构
whos X
size(X)
class(X{1})

% % 显示第一个cell中的第一张图像（假设X{1}包含图像或图像矩阵）
% 加载您的.mat文件（如果尚未加载）
% load('your_file.mat');

% 检查X的基本信息
whos X
size(X)

% 检查X中的第一个元素
if iscell(X)
    % 显示X{1}的类型和大小
    class(X{1})
    size(X{1})
    
    % 如果X{1}也是cell数组，检查更深一层
    if iscell(X{1})
        class(X{1}{1})
        size(X{1}{1})
    end
end

% 检查gt的结构
whos gt
size(gt)
% 显示gt的前几个值以了解其范围和类型
gt(1:min(5, length(gt)))

for i = 1:length(X)
    if isnumeric(X{i}) && ndims(X{i}) >= 3
        % 显示每个cell中的第一张图像
        figure;
        imagesc(X{i}(:,:,1));
        colormap gray;
        title(['First image in X{' num2str(i) '}']);
        axis image;
    end
end