clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=VideoReader('D:\KONGYU\se_2\yiceshi\cs0409.mp4');    %读取视频到avi
t1=zeros(350,1); 
t4=t1;
KUAN_4_CHA1=zeros(350,4);
MAXKUANG_4MJC1=zeros(350,2);
KUAN_4_MJC1=zeros(350,2);
TAR_LOS_TIMES=0;
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


hang_start=350;  %20170329_飞行器，去掉前50帧，通过
hang_end=390;
lie_start=225;                                                                                                                                                                                                                                                                                    00;
lie_end=265;

% hang_start=236;  %20170329_飞行器2，去掉前50帧，通过
% hang_end=253;
% lie_start=448;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;

% 
% hang_start=341;  %20170329_飞行器2cuo，去掉前50帧，通过
% hang_end=383;
% lie_start=414;                                                                                                                                                                                                                                                                                    00;
% lie_end=445;

% hang_start=349;  %20170401_clip，去掉前50帧，通过
% hang_end=375;
% lie_start=420;                                                                                                                                                                                                                                                                                    00;
% lie_end=442;

%灰度拉伸上边界量，值越大则灰度拉伸的范围越窄
HUIDUC=2;
%结束


%追踪区的宽度，根据手动绘制框的宽度来确定
KUANG_2=round(max((hang_end-hang_start),(lie_end-lie_start))*0.5);

% %精确追踪框的宽度
JINGQUE_KUAN=round(KUANG_2*0.5);

K_JINGQUE=KUANG_2-JINGQUE_KUAN;

%膨胀参数
% PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/20);
PENGZHANG=3;
%膨胀结束

PIXEL_NUM=60;  %识别目标的最少像素点个数
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*0.5);%识别目标的最大像素点个数



ERZHIHUA=0.08;  %预设二值化的阈值，这是经验值，在1km目标识别中是有效的，不必更改
ERZHIHUA2=0.8;

YOUXIAO=1;   %目标有效
TAR_LOST=0;  %目标未丢失
TAR_LOST_TIMES=0;   %目标丢失帧数统计
TAR_LOST_SIGN=5;  %确认目标已丢失的帧数上限，当有TAR_LOST_SIGN帧目标均丢失时，认为目标确实已经丢失，请求外部决断

TAR_ON_HOLD=0;

QUEREN=0;

%精确追踪区宽度变量，自适应波门
%%
%分解每一帧存储到pixels，提取行数，列数，帧数
for i=1:avi.NumberOfFrames-50  %循环提取每一帧
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
for i =1 : nFrames-50     
   %%
    if i==1      %判断是否为第一帧图像
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);

%          figure(2);
%          imshow(d(:, :, i));
%          title('差值图1');
%          
       %%
       %划定目标存在区
       dtemp(:,:,i)=d(hang_start:hang_end,lie_start:lie_end,i);
%        figure(3);
%        imshow(dtemp(:, :, i));
%        title('目标区差值图');
       %划定目标区结束
       %%
       %对目标区做灰度值拉伸，输入图像dtemp(:,:,i)，输出图像piclz
       %将图像转换成一列
        picl=dtemp(:);
        SUB_GRAY_AVR_LAST=sum(picl)/((hang_end-hang_start)*(lie_end-lie_start));
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
            d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %求当前帧和下一帧的差，d为求出的差值图像
            d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
            d(:, :, i)=d(:, :, i)+d2(:, :, i);

%             figure(5);
%             imshow(d(:, :, i));
%             title('差值图');



           %提高中心权限  

            %划定中心权限区,输入d，输出bw4
            bw4=d(hangzhong-JINGQUE_KUAN:hangzhong+JINGQUE_KUAN,liezhong-JINGQUE_KUAN :liezhong+JINGQUE_KUAN,i);
            %划定追踪范围结束

            
            figure(6);
            imshow(bw4);
            title('中心权限区');


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
            
            
            figure(7);
            imshow(piclz4);
            title('中心权限区灰度拉伸');
            
            
            picl14=piclz4(:);
            %取出中值
            zhong14=median(picl14)
            if(zhong14<211)
                zhongz14=double(zhong14)*1.2;
            else
                zhongz14=double(zhong14);
            end
            %取出中值结束

            %指数变换，结果为piclz
            shangbian14=zhongz14/255;
            xiabian14=((zhongz14)/255)+((1-zhongz14/255)/HUIDUC);
            if(xiabian14>1)
                xiabian14=1;
            end
            piclz14 = imadjust(piclz4,[shangbian4 xiabian14],[]);

            %提高中心权限结束


            figure(7);
            imshow(piclz14);
            title('中心权限区二次灰度拉伸');

            %划定追踪范围,输入d，输出bw2
            bw2=d(hangzhong-KUANG_2:hangzhong+KUANG_2,liezhong-KUANG_2 :liezhong+KUANG_2,i);
            bw5=bw2;
            
    %         bw7(K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1),K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1))=bw7+piclz4;
            for z=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                for y=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                    bw5(z,y)=piclz14(z-K_JINGQUE+1,y-K_JINGQUE+1);
                end
            end
            bw6=bw5;
            %划定追踪范围结束


            figure(8);
            imshow(bw6);
            title('追踪区中心权限提高完毕');
            
            %s1,按绘图4号框求灰度均值gray_sum;
            tongjiqu=pixel(HUITU_HANG:(HUITU_HANG+HUITU_CHANG),HUITU_LIE:(HUITU_LIE+HUITU_KUAN),i);
            gray_avr=sum(tongjiqu(:))/(HUITU_KUAN*HUITU_CHANG);


          %%
           %对目标区做灰度值拉伸，输入图像bw2，输出图像piclz2
           %将图像转换成一列
            picl=bw2(:);
            picl6=bw6(:);
            %图像转换成一列结束

            %s2,按追踪框求出灰度均值sub_gray_sum
            sub_gray_avr=sum(picl)/((KUANG_2*2+1)*(KUANG_2*2+1));
            %s3,判断追踪区灰度均值和上次的灰度均值差距多少
%xw             
                t=abs(SUB_GRAY_AVR_LAST-sub_gray_avr);
                t1(i,:)=t;
                t2=median(t1);
                t3=max(t1);
                t4(i,:)=t3;
                t5=sum(t4(:))/i;
             
            if ((t>t2)&&(t<t5))%追踪区相减后灰度均值相比
%                if (t<t5)%追踪区相减后灰度均值相比
                %物体在运动
                SUB_GRAY_AVR_LAST=sub_gray_avr;
                YOUXIAO=1;%是否目标区差值图满足条件，如满足（即认为目标确实出现了足够的移动），则执行后续操作
                TAR_ON_HOLD=0;
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
                
                
                figure(9);
                imshow(bw3);
                title('中心权限区二次灰度拉伸后二值化');
                
                
                %膨胀
                se2 = strel('disk',PENGZHANG);
                pengzhang = imdilate(bw3,se2);
                %膨胀结束
                
                figure(10);
                imshow(pengzhang);
                title('中心权限区二次灰度拉伸后二值化后膨胀');

                %求连通域
                [L,n]=bwlabel(pengzhang);   %划定上次追踪，目标位置的KUANG_2*2的方框为追踪区
                STATS =regionprops(L,'Area');
                aa=struct2cell(STATS);
                bb=cell2mat(aa);
                cc=find(bb>PIXEL_NUM&bb<PIXEL_NUM_MAX);
                tf=ismember(L,cc);
                %求连通域，划定追踪范围结束
                 if length(cc)==0
                     YOUXIAO=0;
                 else
                    [L1,n]=bwlabel(tf);
                    STATS1 =regionprops(L1,'Area');
                    aa1=struct2cell(STATS1);
                    bb1=cell2mat(aa1);
                    [max_1 max_loction]=max(bb1);  
                 end

            else%当前可能静止或丢失
                YOUXIAO=0;
                TAR_ON_HOLD=1;
            end
%             if YOUXIAO==0
%                
%             end
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
            hangzhong=hangzhongx+hangzhong-KUANG_2-1;
            liezhong=liezhongy+liezhong-KUANG_2-1;
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
%         if(TAR_ON_HOLD==0)
            if (((GRAY_AVR*1.6)>gray_avr)&&((GRAY_AVR*0.6)<gray_avr))
               TAR_LOST=0;    %目标是否丢失，丢失置1
               GRAY_AVR=gray_avr;
               TAR_ON_HOLD=1;   %如果是这里，即第一次发现目标onhold，则置一，初始状态为0，表示没有静止，若为2，则表明是至少第2次发现目标onhold
               TAR_ON_HOLD_TIMES=1;
              % msgbox('hold');
%                pause;
            else
%                if(TAR_LOST==1)
                   TAR_LOST_TIMES=TAR_LOST_TIMES+1;
                   TAR_LOST=1;
                   msgbox('lost');
%                    pause;
%                end
            end
%         end
    end
     
    if i==1
        if(length(pbb2)~=0)
            HUITU_LIE=lie_start+pbb2(1,1);
            HUITU_HANG=hang_start+pbb2(1,2);
            HUITU_KUAN=pbb2(1,3);
            HUITU_CHANG=pbb2(1,4);
            HUITU_MIANJI=HUITU_HANG*HUITU_KUAN;
            
            tongjiqu=pixel(HUITU_HANG:(HUITU_HANG+HUITU_CHANG),HUITU_LIE:(HUITU_LIE+HUITU_KUAN),i);
            GRAY_AVR=sum(tongjiqu(:))/(HUITU_KUAN*HUITU_CHANG);
        end
    else
        if((length(pbb2)~=0)&&(YOUXIAO==1))
            HUITU_MIANJI=HUITU_CHANG*HUITU_KUAN;
              
            huitu_lie_t=liezhong-KUANG_2+pbb2(1,1);
            huitu_hang_t=hangzhong-KUANG_2+pbb2(1,2);
            huitu_kuan_t=pbb2(1,3);
            huitu_chang_t=pbb2(1,4);
    
           if(abs(huitu_kuan_t*huitu_chang_t-HUITU_MIANJI)<(HUITU_MIANJI*0.8))
             KUAN_4_CHA1(i,:)= [abs(HUITU_KUAN-huitu_kuan_t),abs(HUITU_CHANG-huitu_chang_t)]
             KUAN_4_MJC1(i,:)=[((KUAN_4_CHA1(i,1)+1)*(KUAN_4_CHA1(i,2)+1))]
             KUAN_4_CHA2=[abs(HUITU_LIE-huitu_lie_t), abs(HUITU_HANG-huitu_hang_t)]
             KUAN_4_MJC2(i,:)=[((KUAN_4_CHA1(i,1)+1)*(KUAN_4_CHA1(i,2)+1))]
             MAXKUANG_4MJC1=max(KUAN_4_MJC1(1:i-1))
           if (KUAN_4_MJC1(i,:))<(MAXKUANG_4MJC1(i,:))
               
                HUITU_LIE=huitu_lie_t;
                HUITU_HANG=huitu_hang_t;
                HUITU_KUAN=huitu_kuan_t;
                HUITU_CHANG=huitu_chang_t;
         
              end
        end
        end
    end
     

    
    if (TAR_LOST_TIMES==TAR_LOST_SIGN)
        msgbox('target lost or target stop !')
    end
    figure(1);     %绘图
    imshow(pixel(:, : ,i));
    fprintf('i= %d\n',i);
    rectangle('position',[HUITU_LIE,HUITU_HANG HUITU_KUAN HUITU_CHANG],'edgecolor','r');
    
    
    
    %测试，删
%     figure(2);     %绘图
%     imshow(d(:, : ,i));
%     rectangle('position',[HUITU_LIE HUITU_HANG HUITU_KUAN HUITU_CHANG],'edgecolor','r');
    %测试结束
    
 
%     YOUXIAO=0;
    if i>58
       a=1;                                        
   end
end