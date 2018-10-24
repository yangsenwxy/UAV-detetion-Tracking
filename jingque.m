clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=mmreader('W:\KONGYU\se_2\test_file\cs.mp4');    %读取视频到avi

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
% lie_start=379;                                                                                                                                                                                                                                                                                    00;
% lie_end=397;
% %结束


% hang_start=323;  %20170329_飞行器，去掉前50帧，通过
% hang_end=342;
% lie_start=363;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;


% hang_start=235;  %20170329_飞行器2，去掉前50帧，通过
% hang_end=253;
% lie_start=445;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;

hang_start=253;  %cs，去掉前50帧，通过
hang_end=280;
lie_start=405;                                                                                                                                                                                                                                                                                    00;
lie_end=420;

%灰度拉伸上边界量，值越大则灰度拉伸的范围越窄
HUIDUC=2;
%结束


%追踪区的宽度，根据手动绘制框的宽度来确定
ZHUIZONG_KUAN=round(max((hang_end-hang_start),(lie_end-lie_start))*0.85);

% %精确追踪框的宽度
% JINGQUE_KUAN=round(ZHUIZONG_KUAN*1.25);

%膨胀参数
PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/8);
%膨胀结束

PIXEL_NUM=2;  %识别目标的最少像素点个数
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*1.2);%识别目标的最大像素点个数


ERZHIHUA=0.08;  %预设二值化的阈值，这是经验值，在1km目标识别中是有效的，不必更改

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

%%
%大循环，处理完整个视频才跳出循环
for i =1 : nFrames-100      
   %%
    if i==1      %判断是否为第一帧图像
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);
%        d100(:, :, i) = imsubtract(pixel(:,:,i+1),pixel(:,:,i));   %求当前帧和下一帧的差，d为求出的差值图像
%        d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
%        d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);

         figure(2);
         imshow(d(:, :, i));
         title('差值图1');
         
%          figure(100);
%          imshow(d100(:, :, i));
%          title('差值图100');
%      
       %%
       %划定目标存在区
       dtemp(:,:,i)=d(hang_start:hang_end,lie_start:lie_end,i);
       figure(41);
       imshow(dtemp(:, :, i));
       title('目标区差值图');
       %划定目标区结束
       %%
       %对目标区做灰度值拉伸，输入图像dtemp(:,:,i)，输出图像piclz
       %将图像转换成一列
        picl=dtemp(:);
%         figure(60);
%         imshow(picl);
%         title('一列');
        %图像转换成一列结束

        %取出中值
        zhong=median(picl)
        if(zhong<182)
            zhongz=double(zhong)*1.4;
        else
            zhongz=double(zhong);
        end
        %取出中值结束

        %指数变换，结果为piclz
        shangbian=zhongz/255;
        xiabian=((zhongz)/255)+((1-zhongz/255)/HUIDUC);
        if(xiabian>1)
            xiabian=1;
        end
        piclz = imadjust(dtemp(:,:,i),[shangbian xiabian],[ ]);
%         figure(13);
%         imshow(piclz);
%         title('指数变换后结果图');
        %指数变换结束
       
       
       %对目标区做灰度值拉伸结束
       %%
       
       %无用的语句，只是为了调试时显示
       dtemp2= pixel(hang_start:hang_end,lie_start:lie_end,i);
%        figure(3);
%        imshow(dtemp2);
%        title('划定的目标存在区');   
       %无用的语句到这里结束

       %以最终的二值化阈值再二值化一次
       bw(:, :, i) = im2bw(piclz,ERZHIHUA);              %对差值图像二值化，得到的也就是最终用于处理的图像

       %膨胀
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw(:, :, i),se2);
       %膨胀结束

%        figure(71);
%        imshow(pengzhang);
%        title('膨胀后');

      %以下语句是在求符合条件的区域，将这些区域置1（可能有多个满足条件的区域，最后的图是tf）
      [L,n]=bwlabel(pengzhang);   %求连通域，n为连通域的个数，L为与bw大小相同的矩阵，连通域已标记为（1,2,3....）      
      STATS =regionprops(L,'Area'); %stats是一个长度为L的结构数组（每位是一个数），其中记录的是相应位的面积大小
      aa=struct2cell(STATS);        %结构数组转化为元胞数组（每位是一个数）
      bb=cell2mat(aa);              %将元胞数组转化为矩阵，一行n(连通域个数)列,第一列存储的就是连通域1的面积，第二列连通域2的面积，类推。。
      cc=find((bb>PIXEL_NUM)&(bb<PIXEL_NUM_MAX));         %找出面积满足要求的连通域，返回的是满足要求的面积在bb中的位置,这个位置对应的是其实就是连通域的标号
      %二值化结束

        tf=ismember(L,cc);  %判断L中的元素有没有在cc中出现，tf和index都是和L同样大小的矩阵，
        %tf中为相同位置1，不同置0，这里其实就是将不符合条件的矩阵置0，只保留符合条件的矩阵
% %         figure(4);
% %         imshow(tf);
% %         title('在划定的目标区内二值化的结果');        
       %%
        %对上面处理过的图像tf处理，找出其中面积最大的连通域，认为是目标移动造成的，
        [L1,n]=bwlabel(tf);
        STATS1 =regionprops(L1,'Area');
        aa1=struct2cell(STATS1);
        bb1=cell2mat(aa1);
        [max_1 max_loction]=max(bb1); 
        [tf1 index]=ismember(L1,max_loction);
         
         
%         figure(5);
%         imshow(L1);
%         title('去掉过大过小值后剩余的面积');  
      
    else 
       %%
        d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
        d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
        d(:, :, i)=d(:, :, i)+d2(:, :, i);
        
        
%         figure(2);
%         imshow(d(:, :, i));
%         title('差值图');

        %精确追踪改进区   
        %划定追踪范围,输入d，输出bw2
        bw2=d(hangzhong-ZHUIZONG_KUAN:hangzhong+ZHUIZONG_KUAN,liezhong-ZHUIZONG_KUAN :liezhong+ZHUIZONG_KUAN,i);
        %划定追踪范围结束
       
        
      %%
       %对目标区做灰度值拉伸，输入图像bw2，输出图像piclz2
       %将图像转换成一列
        picl=bw2(:);
%         figure(64);
%         imshow(picl);
%         title('一列');
        %图像转换成一列结束

        %取出中值
        zhong=median(picl)
        if(zhong<211)
            zhongz=double(zhong)*1.2;
        else
            zhongz=double(zhong);
        end
        %取出中值结束

        %指数变换，结果为piclz
        shangbian=zhongz/255;
        xiabian=((zhongz)/255)+((1-zhongz/255)/HUIDUC);
        if(xiabian>1)
            xiabian=1;
        end
        piclz2 = imadjust(bw2,[shangbian xiabian],[]);
%         figure(13);
%         imshow(piclz2);
%         title('指数变换后结果图');
       %指数变换结束
       %对目标区做灰度值拉伸结束
       %%

        bw3 = im2bw(piclz2,ERZHIHUA); 
         
%         figure(19);
%         imshow(bw3);
%         title('二值化化后的区域');
        
       %膨胀
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw3,se2);
       %膨胀结束

%        figure(71);
%        imshow(pengzhang);
%        title('膨胀后');
        
        
        %求连通域
        [L,n]=bwlabel(pengzhang);   %划定上次追踪，目标位置的ZHUIZONG_KUAN*2的方框为追踪区
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>PIXEL_NUM&bb<PIXEL_NUM_MAX);
        tf=ismember(L,cc);
        %求连通域，划定追踪范围结束
        
%         figure(3);
%         imshow(L);
%         title('划定的目标存在区');
%         
%         figure(4);
%         imshow(tf);
%         title('所有联通量');
        
        [L1,n]=bwlabel(tf);
        STATS1 =regionprops(L1,'Area');
        aa1=struct2cell(STATS1);
        bb1=cell2mat(aa1);
       [max_1 max_loction]=max(bb1);  
    end
    %%
%    判断是不是有连通分量，如果有的话，推算连通分量的出真实坐标
    if length(max_loction)   %判断是不是有符合条件的连通分量
         tf1=ismember(L1,max_loction);   %将最大连通分量的位置找出来

%%
        %寻找质心,输入为图像存在唯一区域的tf1
        [L4,n]=bwlabel(tf1);
        STATS =regionprops(L4,'Centroid');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        liezhongy=round(bb(1,1));
        hangzhongx=round(bb(1,2));
        %寻找质心结束，输出为质心（hangzhongx,liezhongy）
%                 
%         figure(20);
%         imshow(L4);
%         title('最大联通分量');
%         text(liezhongy,hangzhongx,'+','color','r');
        %还原真实坐标
        if(i==1)
            hangzhong=hangzhongx+hang_start;
            liezhong=liezhongy+lie_start;
        else
            hangzhong=hangzhongx+hangzhong-ZHUIZONG_KUAN-1;
            liezhong=liezhongy+liezhong-ZHUIZONG_KUAN-1;
        end
        %还原真实坐标结束
        %测试代码，删
%         figure(26);
%         imshow(pixels(:, :, :,i));
%         title('原图');
%         text(liezhong,hangzhong,'+','color','r');
        %测试代码结束
    end
    %%  
 
    [l,m]=bwlabel(tf1);
    status=regionprops(l,'BoundingBox');
    
    
    %测试代码，删
%     figure;
%     imshow(tf1);
%     title('框选测试')
%     rectangle('position',status.BoundingBox,'edgecolor','r');
    %测试代码结束
    
    
    
    paa2=struct2cell(status);
    pbb2=cell2mat(paa2);
    
    
    if i==1
        if(length(pbb2)~=0)
            huitu_lie=lie_start+pbb2(1,1);
            huitu_hang=hang_start+pbb2(1,2);
            huitu_kuan=pbb2(1,3);
            huitu_chang=pbb2(1,4);
        end
    else
        if(length(pbb2)~=0)
            huitu_lie=liezhong-ZHUIZONG_KUAN+pbb2(1,1);
            huitu_hang=hangzhong-ZHUIZONG_KUAN+pbb2(1,2);
            huitu_kuan=pbb2(1,3);
            huitu_chang=pbb2(1,4);
        end
    end
    
    
    figure(1);     %绘图
    imshow(pixel(:, : ,i));
%     hold on
    rectangle('position',[huitu_lie huitu_hang huitu_kuan huitu_chang],'edgecolor','r');
%     text(huitu_lie+huitu_kuan/2,huitu_hang+huitu_chang/2,'+','color','r');
    
%     figure(1);     %绘图
%     imshow(pixels(:, :, :,i));
%     hold on
%     rectangle('Position',[lei_min_last1 hang_min_last1 leikaundu hangkuandu], 'EdgeColor', 'r', 'LineWidth', 2);   %从点lei_min-5 hang_min-8开始画矩形，长为20，宽为20
% %     rectangle('Position',[lei_min hang_min hang_max-hang_min+20 lei_max-lei_min+20], 'EdgeColor', 'r', 'LineWidth', 2);
    a=i;
    
end


% d(:,:,5)=(abs(pixel(:,:,5)-pixel(:,:,1)));
% bw(:,:,5)=im2bw(d(:,:,5),0.13);
% imshow(bw(:,:,5))