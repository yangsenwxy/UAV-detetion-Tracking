clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=mmreader('W:\KONGYU\se_2\test_file\20170401_clip.mp4');    %读取视频到avi

%手动画框的边界,封装后函数的输入变量
% hang_start=274;   %新飞行器_格式_短,去掉前50帧，通过
% hang_end=301;
% lie_start=248;                                                                                                                                                                                                                                                                                    00;
% lie_end=283;

% hang_start=204;%飞行器_固定翼_格式化_1，50，去掉前50帧,通过,
% hang_end=228;
% lie_start=307;                                                                                                                                                                                                                                                                                    00;
% lie_end=336;

% hang_start=308;%1公里人_格式，20，去掉前50帧，通过
% hang_end=339;
% lie_start=337;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;
%结束


% hang_start=323;  %20170329_飞行器，去掉前50帧，通过
% hang_end=342;
% lie_start=363;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;

% hang_start=236;  %20170329_飞行器2，去掉前50帧，通过
% hang_end=253;
% lie_start=448;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;


hang_start=341;  %20170329_飞行器2，去掉前50帧，通过
hang_end=383;
lie_start=414;                                                                                                                                                                                                                                                                                    00;
lie_end=445;

%灰度拉伸上边界量，值越大则灰度拉伸的范围越窄
HUIDUC=2;
%结束


%追踪区的宽度，根据手动绘制框的宽度来确定
ZHUIZONG_KUAN=round(max((hang_end-hang_start),(lie_end-lie_start))*0.75);

% %精确追踪框的宽度
% JINGQUE_KUAN=round(ZHUIZONG_KUAN*1.25);

%膨胀参数
PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/10);
%膨胀结束

PIXEL_NUM=2;  %识别目标的最少像素点个数
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*0.9);%识别目标的最大像素点个数


ERZHIHUA=0.1;  %预设二值化的阈值，这是经验值，在1km目标识别中是有效的，不必更改

%精确追踪区宽度变量，自适应波门
%%
%分解每一帧存储到pixels，提取行数，列数，帧数
for i=1:avi.NumberOfFrames-100  %循环提取每一帧
    img=read(avi,i);           %读取当前帧
    pixels(:,:,:,i)=img;       %将当前帧存储到四维矩阵中的相应层（第i层）去
end

nFrames=size(pixels,4);        %返回矩阵的第四维的size，帧数
rows=size(pixels,1);           %返回矩阵的第一维的size，行
cols=size(pixels,2);           %返回矩阵的第二维的size，列

for i=1:nFrames-50
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i+50)));  %RGB转灰度图像，存入pixel
%   figure(2);
%   imshow(pixels(:,:,:,1));
%   title('原图');
end

k=1;



for i =1 : nFrames-100      
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);
%        d100(:, :, i) = imsubtract(pixel(:,:,i+1),pixel(:,:,i));   %求当前帧和下一帧的差，d为求出的差值图像
%        d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
%        d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);

         figure(2);
         imshow(d(:, :, i));
         title('差值图1');
end