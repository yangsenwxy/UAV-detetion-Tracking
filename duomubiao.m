clear all;          %清除所有变量
clc;                %清除命令窗口
close all           %关闭所有窗口
avi=VideoReader('W:\KONGYU\se_2\test_file\新飞行器_格式_短.mp4');    %通过测试
%avi=VideoReader('W:\KONGYU\se_2\test_file\人体追踪_短.mp4');

%灰度拉伸上边界量，值越大则灰度拉伸的范围越窄
HUIDUC=3;
%结束
MIN_PIX=120;  %设定用于对结果中像素过少的点滤除的阈值，注意这里滤除的是闭操作，膨胀后的图像
TAG_NUM=20;  %设定追踪目标个数
BGN=10;      %用于提取背景的图片数
REAL_TAG_NUM=0;   %对图像处理后，最终分解出的目标个数，可以直接输出
tag_num =TAG_NUM;  %中间参数
ERZHIHUA=0.5;
%%
%分解每一帧存储到pixels，提取行数，列数，帧数
for i=1:avi.NumberOfFrames-100  %循环提取每一帧
    img=read(avi,i);           %读取当前帧
    pixels(:,:,:,i)=img;       %将当前帧存储到四维矩阵中的相应层（第i层）去
end

nFrames=size(pixels,4);        %返回矩阵的第四维的size，帧数
rows=size(pixels,1);           %返回矩阵的第一维的size，行
cols=size(pixels,2);           %返回矩阵的第二维的size，列

for i=1:avi.NumberOfFrames-100
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i)));  %RGB转灰度图像，存入pixel
end

for i=1:avi.NumberOfFrames-15;
    if i==1
        d(:,:)=double(pixel(:,:,i));
    end
    if i<BGN
        d(:,:)= (d(:,:)+ double(pixel(:,:,i)));
        if i==(BGN-1)
            BACG3(:,:)=d(:,:)./(BGN-1);
        end
        i=i+1;
    else
        BACG3(:,:)=round(BACG3(:,:).*39+double(pixel(:,:,i)))/40;
        
        d2(:, :) = imsubtract(uint8((pixel(:,:,i))),uint8(BACG3(:,:)));
        
        figure(19);
        imshow(uint8(BACG3(:,:)));
        title('背景')
        
        figure(7);
        imshow(d2(:,:));
        title('减去背景')
    %     d2(:, :) = imsubtract(pixel(:,:,i+1),BACG3(:,:));
    %     figure(8);
    %     imshow(d2(:,:));
    %     title('减背景2')
    %     d3(:, :) = imsubtract(d2(:, :) , d1(:, :));
    %     figure(9);
    %     imshow(d3(:,:));
    %     title('相减')
    %     d4(:, :) = uint8((double(d2(:, :))+double(d1(:, :)))/2);
    %     figure(10);
    %     imshow(d4(:,:));
    %     title('相加');
    %     d10(:,:)=pixel(:,:,i)-pixel(:,:,i+1);
    %     figure(17);
    %     imshow(d10(:,:));
    % %     title('相加');
    %     d4(:,:)=d4(:, :)+8*d10(:,:);
    %      figure(18);
    %     imshow(d4(:,:));
    %     title('相加');
        %%
        %滤波
        d5=filter2(fspecial('average',4),d2)/255;          %进行5*5模板平滑滤波
        figure(11);
        imshow(d5);
        title('滤波后')
    %%
    %灰度拉伸
     %对目标区做灰度值拉伸，输入图像d5，输出图像piclz
           %将图像转换成一列
            picl=d5(:);
    %         figure(60);
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
            piclz = imadjust(d5,[shangbian xiabian],[]);
            figure(13);
            imshow(piclz);
            title('指数变换后结果图');
            %指数变换结束   
    %灰度拉伸结束

     %二值化可信度
       if i<20
           ERZHIHUA=0.5;
       else if i<40
                ERZHIHUA=0.4;
           else if i<60
                    ERZHIHUA=0.3;
               else
                    ERZHIHUA=0.2;
               end
           end
       end
       %判断可信度结束
        %二值化
        d6 = im2bw(piclz, ERZHIHUA);
        figure(12);
        imshow(d6);
        title('二值化')
        %%
        %边缘提取
    %     bianyuantu=edge(d6,'prewitt');%用prewitt算子进行边缘检测
    %     figure;
    %     imshow(bianyuan);
    %     title('边缘')
    %     d6(1:100,:)=0;
        se=[0,0,1,0,0;0,1,1,1,0;1,1,1,1,1;0,1,1,1,0;0,0,1,0,0];
        pengzhang = imclose(d6,se); %闭操作
        se2 = strel('disk',3);
        pengzhang2 = imdilate(pengzhang,se2);  %膨胀
    %     
    %     pengzhang(1,:)=0;
    %     pengzhang(rows,:)=0;
    %     pengzhang(:,1)=0;
    %     pengzhang(:,cols)=0;
        figure(13);
        imshow(pengzhang2);
        title('膨胀后');

       %%
       %像素点可信度
       if i<20
           MIN_PIX=120;
       else if i<40
               MIN_PIX=100;
           else if i<60
                   MIN_PIX=80;
               else
                   MIN_PIX=60;
               end
           end
       end
       %判断可信度结束

       %%
        %求连通域
         figure(5);
         imshow(pixel(:,:,i));
         [L,n]=bwlabel(pengzhang2);   %求连通域，n为连通域的个数，L为与bw大小相同的矩阵，连通域已标记为（1,2,3....）

    %          figure(22);
    %          imshow(L)
    %          title('所有连通域')
            if(n~=0)
                STATS =regionprops(L,'Area'); %stats是一个长度为L的结构数组（每位是一个数），其中记录的是相应位的面积大小
                aa=struct2cell(STATS);        %结构数组转化为元胞数组（每位是一个数）
                bb=cell2mat(aa);              %将元胞数组转化为矩阵，一行n(连通域个数)列,第一列存储的就是连通域1的面积，第二列连通域2的面积，类推。。
                cc=find(bb>MIN_PIX);         %找出面积满足要求的连通域，返回的是满足要求的面积在bb中的位置,这个位置对应的是其实就是连通域的标号
        %          if length(cc)
                if(cc~=0)
                    tf=ismember(L,cc);  %判断L中的元素有没有在cc中出现，tf是和L同样大小的矩阵，
            %           %tf中为相同位置1，不同置0，这里其实就是将不符合条件的矩阵置0，只保留符合条件的矩阵
    %                  figure(21);
    %                  imshow(tf);
    %                  title('去除过大过小区域的连通域');
                     %到这里判断是否是符合条件的区域结束

                     %%
                     %取出符合预计的n（tag_num）个目标区域
                    [L1,n]=bwlabel(tf);
                    if(n~=0)
                        STATS1 =regionprops(L1,'Area'); 
                        aa1=struct2cell(STATS1);
                        bb1=cell2mat(aa1);
                        num_can_do=length(bb1);%多目标修改区1
                        if(num_can_do<tag_num)
                            tag_num=num_can_do;
                        end
                        for tag_num_i=1:tag_num     
                            [max_1(tag_num_i) max_loction]=max(bb1);%max_1是最大值，max_location是最大值第一次出现的位置,其实就是面积最大连通域的标号
                            tf1(:,:,tag_num_i)=ismember(L1,max_loction);%将面积最大的连通域置1，其他位置置0，得图像tf1，tf1中就只有，存在目标的区域
                            bb1(max_loction)=[1];

                            figure(20);
                            imshow(tf1(:,:,tag_num_i));
                            title('分解出的目标');
                        end
                        tf2=tf1(:,:,1);
                        for tag_num_i=2:tag_num
                            tf2=tf2+tf1(:,:,tag_num_i);
                        end
                        tag_num=TAG_NUM;

                        figure;
                        imshow(tf2);
                        title('所有目标')
                       %%

                         %绘制目标框
                         [l,m]=bwlabel(tf2,8);
                         status=regionprops(l,'BoundingBox');
                         tar_position_temp=struct2cell(status);
                         tar_position=cell2mat(tar_position_temp);    %(列，行，列宽，行宽)

                         REAL_TAG_NUM=(length(tar_position))/4;

                         figure(71);
                         imshow(pixel(:,:,i));
                         title('检测结果');
                         hold on;

                         jishu=1;
                         for tar=1:REAL_TAG_NUM
                             reall=ceil(tar_position(jishu));
                             jishu=jishu+1;
                             realh=ceil(tar_position(jishu));
                             jishu=jishu+1;
                             reallk=tar_position(jishu);
                             jishu=jishu+1;
                             realhk=tar_position(jishu);
                             jishu=jishu+1;
                             rectangle('position',[reall realh reallk realhk],'edgecolor','r');
                         end
                    end
                end
            end
    end
end

