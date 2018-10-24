clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=mmreader('W:\KONGYU\se_2\test_file\20170401_clip.mp4');    %读取视频到avi

%手动画框的边界,封装后函数的输入变量
% hang_start=274;   %新飞行器_格式_短,去掉前50帧，通过
% hang_end=301;
% lie_start=248;                                                                                                                                                                                                                                                                                    00;
% % lie_end=283;
% 
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

% 
% hang_start=341;  %20170329_飞行器2，去掉前50帧，通过
% hang_end=383;
% lie_start=414;                                                                                                                                                                                                                                                                                    00;
% lie_end=445;

hang_start=349;  %20170401_clip，去掉前50帧，通过
hang_end=375;
lie_start=420;                                                                                                                                                                                                                                                                                    00;
lie_end=442;

%灰度拉伸上边界量，值越大则灰度拉伸的范围越窄
HUIDUC=2;
%结束


%追踪区的宽度，根据手动绘制框的宽度来确定
ZHUIZONG_KUAN=round(max((hang_end-hang_start),(lie_end-lie_start))*0.75);

% %精确追踪框的宽度
JINGQUE_KUAN=round(ZHUIZONG_KUAN*0.5);

K_JINGQUE=ZHUIZONG_KUAN-JINGQUE_KUAN;

%膨胀参数
PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/10);
%膨胀结束

PIXEL_NUM=2;  %识别目标的最少像素点个数
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*0.9);%识别目标的最大像素点个数



ERZHIHUA=0.08;  %预设二值化的阈值，这是经验值，在1km目标识别中是有效的，不必更改
ERZHIHUA2=0.8;

YOUXIAO=1;   %目标有效
TAR_LOST=0;  %目标未丢失
TAR_LOST_TIMES=0;   %目标丢失帧数统计
TAR_LOST_SIGN=5;  %确认目标已丢失的帧数上限，当有TAR_LOST_SIGN帧目标均丢失时，认为目标确实已经丢失，请求外部决断

TAR_ON_HOLD=0;

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

         figure(2);
         imshow(d(:, :, i));
         title('差值图1');
         
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
        SUB_GRAY_SUM=sum(picl)/((hang_end-hang_start)*(lie_end-lie_start));
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
        piclz = imadjust(dtemp(:,:,i),[shangbian xiabian],[]);
        %指数变换结束
       %对目标区做灰度值拉伸结束
       %%
       %以最终的二值化阈值再二值化一次
       bw(:, :, i) = im2bw(piclz,ERZHIHUA);              %对差值图像二值化，得到的也就是最终用于处理的图像

       %膨胀
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw(:, :, i),se2);
       %膨胀结束

      %以下语句是在求符合条件的区域，将这些区域置1（可能有多个满足条件的区域，最后的图是tf）
      [L,n]=bwlabel(pengzhang);   %求连通域，n为连通域的个数，L为与bw大小相同的矩阵，连通域已标记为（1,2,3....）      
      STATS =regionprops(L,'Area'); %stats是一个长度为L的结构数组（每位是一个数），其中记录的是相应位的面积大小
      aa=struct2cell(STATS);        %结构数组转化为元胞数组（每位是一个数）
      bb=cell2mat(aa);              %将元胞数组转化为矩阵，一行n(连通域个数)列,第一列存储的就是连通域1的面积，第二列连通域2的面积，类推。。
      cc=find((bb>PIXEL_NUM)&(bb<PIXEL_NUM_MAX));         %找出面积满足要求的连通域，返回的是满足要求的面积在bb中的位置,这个位置对应的是其实就是连通域的标号
      %二值化结束

        tf=ismember(L,cc);  %判断L中的元素有没有在cc中出现，tf和index都是和L同样大小的矩阵，
        %tf中为相同位置1，不同置0，这里其实就是将不符合条件的矩阵置0，只保留符合条件的矩阵
;        
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
        YOUXIAO=1;
      
    else 
       %%
       if (TAR_ON_HOLD~=2)&&(TAR_ON_HOLD~=1)   %目标既未进入静止状态
            d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
            d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
            d(:, :, i)=d(:, :, i)+d2(:, :, i);

    %         figure(2);
    %         imshow(d(:, :, i));
    %         title('差值图');



           %提高中心权限  

            %划定中心权限区,输入d，输出bw2
            bw4=d(hangzhong-JINGQUE_KUAN:hangzhong+JINGQUE_KUAN,liezhong-JINGQUE_KUAN :liezhong+JINGQUE_KUAN,i);
            %划定追踪范围结束


           %对目标区做灰度值拉伸，输入图像bw2，输出图像piclz2
           %将图像转换成一列
            picl4=bw4(:);
            %取出中值
            zhong4=median(picl4)
            if(zhong4<211)
                zhongz4=double(zhong4)*1.2;
            else
                zhongz4=double(zhong4);
            end
            %取出中值结束

            %指数变换，结果为piclz
            shangbian4=zhongz4/255;
            xiabian4=((zhongz4)/255)+((1-zhongz4/255)/HUIDUC);
            if(xiabian4>1)
                xiabian4=1;
            end
            piclz4 = imadjust(bw4,[shangbian4 xiabian4],[]);

            %提高中心权限结束



            %划定追踪范围,输入d，输出bw2
            bw2=d(hangzhong-ZHUIZONG_KUAN:hangzhong+ZHUIZONG_KUAN,liezhong-ZHUIZONG_KUAN :liezhong+ZHUIZONG_KUAN,i);
            bw5=bw2;
    %         bw7(K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1),K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1))=bw7+piclz4;
            for z=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                for y=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                    bw5(z,y)=piclz4(z-K_JINGQUE+1,y-K_JINGQUE+1);
                end
            end
            bw6=bw5;
            %划定追踪范围结束


            %s1,按上次小框求灰度均值gray_sum;
            tongjiqu=pixel(HUITU_HANG:(HUITU_HANG+HUITU_CHANG),HUITU_LIE:(HUITU_LIE+HUITU_KUAN),i);
            gray_sum=sum(tongjiqu(:))/(HUITU_KUAN*HUITU_CHANG);


          %%
           %对目标区做灰度值拉伸，输入图像bw2，输出图像piclz2
           %将图像转换成一列
            picl=bw2(:);
            picl6=bw6(:);
            %图像转换成一列结束

            %s2,按追踪框求出灰度均值sub_gray_sum
            sub_gray_sum=sum(picl)/((ZHUIZONG_KUAN+1)*(ZHUIZONG_KUAN+1));
            %s3,判断追踪区灰度均值和上次的灰度均值差距多少
            if abs(SUB_GRAY_SUM-sub_gray_sum)<3  %追踪区相减后灰度均值相比
                SUB_GRAY_SUM=sub_gray_sum;
                YOUXIAO=1;%是否目标区差值图满足条件，如满足（即认为目标确实出现了足够的移动），则执行后续操作
                %取出中值
                zhong=median(picl6);
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
                piclz2 = imadjust(bw6,[shangbian xiabian],[]);
               %指数变换结束
               %对目标区做灰度值拉伸结束
%%

                bw3 = im2bw(piclz2,ERZHIHUA); 
                
                %膨胀
                se2 = strel('disk',PENGZHANG);
                pengzhang = imdilate(bw3,se2);
                %膨胀结束

                %求连通域
                [L,n]=bwlabel(pengzhang);   %划定上次追踪，目标位置的ZHUIZONG_KUAN*2的方框为追踪区
                STATS =regionprops(L,'Area');
                aa=struct2cell(STATS);
                bb=cell2mat(aa);
                cc=find(bb>PIXEL_NUM&bb<PIXEL_NUM_MAX);
                tf=ismember(L,cc);
                %求连通域，划定追踪范围结束

                [L1,n]=bwlabel(tf);
                STATS1 =regionprops(L1,'Area');
                aa1=struct2cell(STATS1);
                bb1=cell2mat(aa1);
                [max_1 max_loction]=max(bb1);  

            else
                YOUXIAO=0
            end
       else
           %TAR_ON_HOLD==1
           %执行语句
           if (TAR_ON_HOLD_TIMES==1)
               TAR_ON_HOLD_DOMAIN=pixel(HUITU_HANG-10:(HUITU_HANG+HUITU_CHANG+10),HUITU_LIE-10:(HUITU_LIE+HUITU_KUAN)+10,i);
               TAR_ON_HOLD_EDGE_LAST=im2bw(TAR_ON_HOLD_DOMAIN,ERZHIHUA2);
               TAR_ON_HOLD_TIMES=TAR_ON_HOLD_TIMES+1;
           else
               TAR_ON_HOLD_DOMAIN=pixel(HUITU_HANG-10:(HUITU_HANG+HUITU_CHANG+10),HUITU_LIE-10:(HUITU_LIE+HUITU_KUAN)+10,i);
               TAR_ON_HOLD_EDGE=im2bw(TAR_ON_HOLD_DOMAIN,ERZHIHUA2);
               TAR_ON_HOLD_LOST_DIF=abs(TAR_ON_HOLD_EDGE-TAR_ON_HOLD_EDGE_LAST);
               [tara tarb]=find(TAR_ON_HOLD_LOST_DIF==1);
               tarc=length(tara);
               if tarc<3   %小于25个坐标点，则认为目标没有移动
                  TAR_ON_HOLD_DOMAIN_LAST=TAR_ON_HOLD_DOMAIN;
                  TAR_ON_HOLD=2;
                  YOUXIAO=0;
                  msgbox('TARGET HOLD');
                  pause;
               else
                   TAR_ON_HOLD_TIMES=0;
                   TAR_ON_HOLD=0;
                   YOUXIAO=0;
               end
           end
       end
    end
    %%
%    判断是不是有连通分量，如果有的话，推算连通分量的出真实坐标
    if (YOUXIAO==1)&&(length(max_loction))   %判断是不是有符合条件的连通分量
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
        
        %还原真实坐标
        if(i==1)
            hangzhong=hangzhongx+hang_start;
            liezhong=liezhongy+lie_start;
        else
            hangzhong=hangzhongx+hangzhong-ZHUIZONG_KUAN-1;
            liezhong=liezhongy+liezhong-ZHUIZONG_KUAN-1;
        end
        %还原真实坐标结束

    end
    %% 
    if(YOUXIAO==1)
        [l,m]=bwlabel(tf1);
        status=regionprops(l,'BoundingBox');

        paa2=struct2cell(status);
        pbb2=cell2mat(paa2);
    else
        if(TAR_ON_HOLD==0)
            if (((GRAY_SUM*2)>gray_sum)&&((GRAY_SUM*0.2)<gray_sum))
               TAR_LOST=0;    %目标是否丢失，丢失置1
               GRAY_SUM=gray_sum;
               TAR_ON_HOLD=1;   %如果是这里，即第一次发现目标onhold，则置一，初始状态为0，表示没有静止，若为2，则表明是至少第2次发现目标onhold
               TAR_ON_HOLD_TIMES=1;
            else
               if(TAR_LOST==1)
                   TAR_LOST_TIMES=TAR_LOST_TIMES+1;
                   TAR_LOST=1;
                   GRAY_SUM=GRAY_SUM;
               else
                   TAR_LOST=1;
                   GRAY_SUM=GRAY_SUM;
               end
            end
        end
    end
     
    if i==1
        if(length(pbb2)~=0)
            HUITU_LIE=lie_start+pbb2(1,1);
            HUITU_HANG=hang_start+pbb2(1,2);
            HUITU_KUAN=pbb2(1,3);
            HUITU_CHANG=pbb2(1,4);
            HUITU_MIANJI=HUITU_HANG*HUITU_KUAN;
            
            tongjiqu=pixel(HUITU_HANG:(HUITU_HANG+HUITU_CHANG),HUITU_LIE:(HUITU_LIE+HUITU_KUAN),i);
            GRAY_SUM=sum(tongjiqu(:))/(HUITU_KUAN*HUITU_CHANG);
        end
    else
        if((length(pbb2)~=0)&&(YOUXIAO==1))
            
            HUITU_MIANJI=HUITU_CHANG*HUITU_KUAN;
            
            huitu_lie_t=liezhong-ZHUIZONG_KUAN+pbb2(1,1);
            huitu_hang_t=hangzhong-ZHUIZONG_KUAN+pbb2(1,2);
            huitu_kuan_t=pbb2(1,3);
            huitu_chang_t=pbb2(1,4);
            
            if(abs(huitu_kuan_t*huitu_chang_t-HUITU_MIANJI)<(HUITU_MIANJI*0.3))
                HUITU_LIE=huitu_lie_t;
                HUITU_HANG=huitu_hang_t;
                HUITU_KUAN=huitu_kuan_t;
                HUITU_CHANG=huitu_chang_t;
            end
        end
        
    end
    
    if (TAR_LOST_TIMES==TAR_LOST_SIGN)
        msgbox('target lost!')
        pause;
    end
    figure(1);     %绘图
    imshow(pixel(:, : ,i));
    rectangle('position',[HUITU_LIE HUITU_HANG HUITU_KUAN HUITU_CHANG],'edgecolor','r');

    YOUXIAO=0;
end