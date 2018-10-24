clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=VideoReader('W:\KONGYU\se_2\test_file\飞行器_固定翼_格式化_1.mp4');    %读取视频到avi

ERZHIHUA=0.05;
JINGQUE_L=120;
JINGQUE_2L=JINGQUE_L*2+1;
%%
%分解每一帧存储到pixels，提取行数，列数，帧数
for i=1:avi.NumberOfFrames-200  %循环提取每一帧
    img=read(avi,i);           %读取当前帧
    pixels(:,:,:,i)=img;       %将当前帧存储到四维矩阵中的相应层（第i层）去
end

nFrames=size(pixels,4);        %返回矩阵的第四维的size，帧数
rows=size(pixels,1);           %返回矩阵的第一维的size，行
cols=size(pixels,2);           %返回矩阵的第二维的size，列

for i=1:nFrames
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i)));  %RGB转灰度图像，存入pixel
%   figure;
%   imshow(pixels(:,:,:,1));
%   title('yuantu');
%   figure;
%   imshow(pixel(:,:,i));
%   title('huidu')
end

k=1;

%%
%大循环，处理完整个视频才跳出循环
for i = 1 : nFrames-4       
    d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
    d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
    d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);
    bw(:, :, i) = im2bw(d(:, :, i), ERZHIHUA);                %对差值图像二值化，得到的也就是最终用于处理的图像
%     
%      figure;
%      imshow(d(:,:,i));
%      title('灰度差值图1')
     
%      figure;
%      imshow(d1(:,:,i));
%      title('灰度差值图2')
%      
%      figure;
%      imshow(d3(:,:,i));
%      title('差值图相加')
%      
% %      imhist(d(:,:,i));
%      figure;
%      imshow(pixel(:,:,i))
%      title('灰度图1');
%      figure;
%      imshow(pixel(:,:,i+1))
%      title('灰度图2')
%      figure;
%      imshow(bw(:, :, i));
%      title('二值化的差值图');
% %     
   %%
    if i==1      %判断是否为第一帧图像
        
        
        %以下语句是在求符合条件的区域，将这些区域置1（可能有多个满足条件的区域，最后的图是tf）
        [L,n]=bwlabel(bw(:, :, i));   %求连通域，n为连通域的个数，L为与bw大小相同的矩阵，连通域已标记为（1,2,3....）
        
%          figure;
%          imshow(L)
%          title('所有连通域')
%         
        STATS =regionprops(L,'Area'); %stats是一个长度为L的结构数组（每位是一个数），其中记录的是相应位的面积大小
        aa=struct2cell(STATS);        %结构数组转化为元胞数组（每位是一个数）
        bb=cell2mat(aa);              %将元胞数组转化为矩阵，一行n(连通域个数)列,第一列存储的就是连通域1的面积，第二列连通域2的面积，类推。。
        cc=find(bb>2&bb<100);         %找出面积满足要求的连通域，返回的是满足要求的面积在bb中的位置,这个位置对应的是其实就是连通域的标号
%          if length(cc)
         tf=ismember(L,cc);  %判断L中的元素有没有在cc中出现，tf和index都是和L同样大小的矩阵，
          %tf中为相同位置1，不同置0，这里其实就是将不符合条件的矩阵置0，只保留符合条件的矩阵
%          figure;
%          imshow(tf);
%          title('去除过大过小区域的二值图');
         %到这里判断是否是符合条件的区域结束
         
         
       %%
        %对上面处理过的图像tf处理，找出其中面积最大的连通域，认为是目标移动造成的，
        [L,n]=bwlabel(bw(:, :, i));
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%          if length(cc)
          [tf index]=ismember(L,cc);
%          figure,
%          imshow(tf);
%          title('符合预期目标')
         
        se=[0 0 1 0 0,0 0 1 0 0,1 1 1 1 1,0 0 1 0 0,0 0 1 0 0];    %结构单元
        pengzhang = imclose(tf,se);     %用se
        se2 = strel('disk',10);
        tf= imdilate(tf,se2);
%         figure;
%         imshow(tf);
%         title('闭操作');
         
         
        [L1,n]=bwlabel(tf);
         STATS1 =regionprops(L1,'Area');
         aa1=struct2cell(STATS1);
         bb1=cell2mat(aa1);
         [max_1 max_loction]=max(bb1); 
         [tf1 index]=ismember(L1,max_loction);
        
      
    else 
       %%
        %求连通域，标定大致可能的区域，同上
        [L,n]=bwlabel(bw(:, :, i));
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%       if length(cc)
        tf=ismember(L,cc);
       
        if(((hang_min-20)<1)||((lei_min-20)<1)||((hang_min+20)>rows)||((lei_min+20)>cols))
            msgbox('飞出范围');
        end
        [L,n]=bwlabel(tf(hang_min-20:hang_min+20,lei_min-20 :lei_min+20));   %hang_min和lei_min是上次绘点的参考点，是真实坐标，这时候在上次参考点为中心41*41的范围内寻找目标
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%       if length(cc)
        tf=ismember(L,cc);
%         figure,
%         imshow(tf);
%         title('符合预期目标');
        
       se=[0 0 1 0 0,0 0 1 0 0,1 1 1 1 1,0 0 1 0 0,0 0 1 0 0];    %结构单元
        pengzhang = imclose(tf,se);     %用se
        se2 = strel('disk',10);
        tf= imdilate(tf,se2);
%         figure;
%         imshow(tf);
%         title('闭操作');
        
%       确定是否存在符合预期的目标
       [L1,n]=bwlabel(tf);
       STATS1 =regionprops(L1,'Area');
       aa1=struct2cell(STATS1);
       bb1=cell2mat(aa1);
      [max_1 max_loction]=max(bb1);  
    end
    %%
    %判断是不是有连通分量，如果有的话，推算连通分量的出真实坐标
    if length(max_loction)   %判断是不是有符合条件的连通分量
         tf1=ismember(L1,max_loction);   %将最大连通分量的位置找出来

%%
          %取最左上的点为参考点
          if i==1
              gg=find(tf1==1)';   %找出满足条件的位置标号，存放在一维矩阵gg中，这里gg先行后列计数，大小可能达到cow*lei
              lei=fix(gg./rows)+1; %朝0取整
              lei_min=min(lei);   %取出最小列值
%               lei_min_1=lei_min;
              lei_max=max(lei);
              hang=rem(gg,rows)+1; %取余，最小的行值
              hang_min=min(hang);
%               hang_min_1=hang_min;
              hang_max=max(hang);
          else 
              gg=find(tf1==1)';
              lei=fix(gg./41)+1-21+lei_min_last; %求真实参考坐标
%               lei_min_1=min(lei);
              lei_min=lei(1);
              lei_max=max(lei);
              hang=rem(gg,41)+1-21+hang_min_last;
%               hang_min_1=min(hang);
              hang_min=hang(1);
              hang_max=max(hang);
          end
          
          if k==1
              lei_min_last=lei_min;
              hang_min_last=hang_min;
              k=k+1;
          end
    
%           if ((abs(lei_min-lei_min_last)<100)||(abs(hang_min-hang_min_last)<100))
              lei_min_last=lei_min;
              lei_max_last= lei_max;
              hang_min_last=hang_min;
              hang_max_last=hang_max;
%           else 
%               lei_min=lei_min_last;
%               hang_min=hang_min_last;
%               lei_max=lei_max_last; 
%               hang_max=hang_max_last;
%           end
    end
    %%
    %精确追踪改进区
    %mubiaoqu=d3(hang_min_last-100:hang_min_last+100,lei_min_last-100:lei_min_last+100,i);
     if(((hang_min_last-JINGQUE_L)<1)||((lei_min_last-JINGQUE_L)<1)||((hang_min_last+JINGQUE_L)>rows)||((lei_min_last+JINGQUE_L)>cols))
            msgbox('飞出范围');
     end
      edge_top=hang_min_last-JINGQUE_L;
      edge_bottom=hang_min_last+JINGQUE_L;
      edge_left=lei_min_last-JINGQUE_L;
      edge_right=lei_min_last+JINGQUE_L;
    %%
    %开始分离目标
    mubiaoqu=bw(edge_top:edge_bottom,edge_left:edge_right,i);
%     figure;
%     imshow(mubiaoqu);
%     title('目标区')
    bianyuantu=edge(mubiaoqu,'prewitt');%用prewitt算子进行边缘检测
    
%     figure(3);
%     imshow(bianyuantu);
%     title('边缘')
%     
    se=[0 0 0 0 1 0 0 0 0,0 0 0 0 1 0 0 0 0,0 0 0 1 1 1 0 0 0,0 0 1 1 1 1 1 0 0,1 1 1 1 1 1 1 1 1,0 0 1 1 1 1 1 0 0,0 0 0 1 1 1 0 0 0,0 0 0 0 1 0 0 0 0,0 0 0 0 1 0 0 0 0];    %结构单元
    bicaozuo = imclose(bianyuantu,se);     %用se
    
%     figure;
%     imshow(bicaozuo);
%     title('闭操作')
    
    se2 = strel('disk',20);
    pengzhang = imdilate(bicaozuo,se2);  %膨胀
%     figure(2);
%     imshow(pengzhang);
%     title('膨胀后');
%     
% 
    pengzhang(1,:)=0;
    pengzhang(JINGQUE_2L,:)=0;
    pengzhang(:,1)=0;
    pengzhang(:,JINGQUE_2L)=0;

    %面积符合区
    [PL1,n]=bwlabel(pengzhang);
    STATS2 =regionprops(PL1,'Area'); 
    paa1=struct2cell(STATS2);
    pbb1=cell2mat(paa1);
    cc=find(pbb1>200&pbb1<10000);
    tf3=ismember(PL1,cc);
%     figure(3);
%     imshow(tf3);
%     title('面积区');
%     %最大面积区
    [L4,n]=bwlabel(tf3);
    STATS4 =regionprops(L4,'Area');
    aa4=struct2cell(STATS4);
    bb4=cell2mat(aa4);
    [max_val max_loc1]=max(bb4);
    tf4=ismember(L4,max_loc1);
%     figure(4);
%     imshow(tf4);
%     title('最大区');
    %求出实际坐标值
    weiz=find(tf4==1)';   %找出满足条件的位置标号，存放在一维矩阵gg中，这里gg先行后列计数，大小可能达到cow*lei
    lei=fix(weiz./JINGQUE_2L)+1; %朝0取整
    lei_min1=min(lei)-JINGQUE_L+lei_min;   %取出最小列值
    lei_max1=max(lei)+lei_max-JINGQUE_L;
    hang=rem(weiz,JINGQUE_2L)+1; %取余，最小的行值
    hang_min1=min(hang)+hang_min-JINGQUE_L;
    hang_max1=max(hang)+hang_max-JINGQUE_L;
    leikuandu=max(lei)-min(lei);
    hangkuandu=max(hang)-min(hang);
   %% 
    %确定此次是否找到目标，如果没有找到，则保留上次位置
    if((length(lei_min1)==0)||(length(lei_max1)==0)||(length(hang_min1)==0)||(length(hang_max1)==0))
        lei_min1=lei_min_last1;
%         lei_max_last1= lei_max1;
        hang_min1=hang_min_last1;
%         hang_max_last1=hang_max1;
        leikuandu=leikuandu_last;
        hangkuandu=hangkuandu_last;
    end 
    lei_min_last1=lei_min1;
    hang_min_last1=hang_min1;
    leikuandu_last=leikuandu;
    hangkuandu_last=hangkuandu;
    
%     figure,
%     imshow(tf);
    figure(1);     %绘图
    imshow(pixels(:, :, :, i));
    hold on
    %%
    %标定图像
    rectangle('Position',[lei_min1 hang_min1 leikuandu hangkuandu], 'EdgeColor', 'r', 'LineWidth', 2);   %从点lei_min-5 hang_min-8开始画矩形，长为20，宽为20
%     rectangle('Position',[lei_min hang_min hang_max-hang_min+20 lei_max-lei_min+20], 'EdgeColor', 'r', 'LineWidth', 2);
    a1212a=1;
    
end


% d(:,:,5)=(abs(pixel(:,:,5)-pixel(:,:,1)));
% bw(:,:,5)=im2bw(d(:,:,5),0.13);
% imshow(bw(:,:,5))