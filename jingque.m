clear all;          %������б���
clc;                %��������
close all           %�ر����д���
avi=mmreader('W:\KONGYU\se_2\test_file\cs.mp4');    %��ȡ��Ƶ��avi

%�ֶ�����ı߽�,��װ�������������
% hang_start=274;   %�·�����_��ʽ_��,ȥ��ǰ50֡��ͨ��
% hang_end=301;
% lie_start=248;                                                                                                                                                                                                                                                                                    00;
% lie_end=283;

% hang_start=204;%������_�̶���_��ʽ��_1��50��ȥ��ǰ50֡,ͨ��,
% hang_end=228;
% lie_start=307;                                                                                                                                                                                                                                                                                    00;
% lie_end=336;

% hang_start=308;%1������_��ʽ��20��ȥ��ǰ50֡��ͨ��
% hang_end=339;
% lie_start=379;                                                                                                                                                                                                                                                                                    00;
% lie_end=397;
% %����


% hang_start=323;  %20170329_��������ȥ��ǰ50֡��ͨ��
% hang_end=342;
% lie_start=363;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;


% hang_start=235;  %20170329_������2��ȥ��ǰ50֡��ͨ��
% hang_end=253;
% lie_start=445;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;

hang_start=253;  %cs��ȥ��ǰ50֡��ͨ��
hang_end=280;
lie_start=405;                                                                                                                                                                                                                                                                                    00;
lie_end=420;

%�Ҷ������ϱ߽�����ֵԽ����Ҷ�����ķ�ΧԽխ
HUIDUC=2;
%����


%׷�����Ŀ�ȣ������ֶ����ƿ�Ŀ����ȷ��
ZHUIZONG_KUAN=round(max((hang_end-hang_start),(lie_end-lie_start))*0.85);

% %��ȷ׷�ٿ�Ŀ��
% JINGQUE_KUAN=round(ZHUIZONG_KUAN*1.25);

%���Ͳ���
PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/8);
%���ͽ���

PIXEL_NUM=2;  %ʶ��Ŀ����������ص����
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*1.2);%ʶ��Ŀ���������ص����


ERZHIHUA=0.08;  %Ԥ���ֵ������ֵ�����Ǿ���ֵ����1kmĿ��ʶ��������Ч�ģ����ظ���

%��ȷ׷������ȱ���������Ӧ����
%%
%�ֽ�ÿһ֡�洢��pixels����ȡ������������֡��
for i=1:avi.NumberOfFrames-100  %ѭ����ȡÿһ֡
    img=read(avi,i);           %��ȡ��ǰ֡
    pixels(:,:,:,i)=img;       %����ǰ֡�洢����ά�����е���Ӧ�㣨��i�㣩ȥ
end

nFrames=size(pixels,4);        %���ؾ���ĵ���ά��size��֡��
rows=size(pixels,1);           %���ؾ���ĵ�һά��size����
cols=size(pixels,2);           %���ؾ���ĵڶ�ά��size����

for i=1:nFrames-50
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i+50)));  %RGBת�Ҷ�ͼ�񣬴���pixel
%   figure(2);
%   imshow(pixels(:,:,:,1));
%   title('ԭͼ');
end

k=1;

%%
%��ѭ����������������Ƶ������ѭ��
for i =1 : nFrames-100      
   %%
    if i==1      %�ж��Ƿ�Ϊ��һ֡ͼ��
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);
%        d100(:, :, i) = imsubtract(pixel(:,:,i+1),pixel(:,:,i));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
%        d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
%        d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);

         figure(2);
         imshow(d(:, :, i));
         title('��ֵͼ1');
         
%          figure(100);
%          imshow(d100(:, :, i));
%          title('��ֵͼ100');
%      
       %%
       %����Ŀ�������
       dtemp(:,:,i)=d(hang_start:hang_end,lie_start:lie_end,i);
       figure(41);
       imshow(dtemp(:, :, i));
       title('Ŀ������ֵͼ');
       %����Ŀ��������
       %%
       %��Ŀ�������Ҷ�ֵ���죬����ͼ��dtemp(:,:,i)�����ͼ��piclz
       %��ͼ��ת����һ��
        picl=dtemp(:);
%         figure(60);
%         imshow(picl);
%         title('һ��');
        %ͼ��ת����һ�н���

        %ȡ����ֵ
        zhong=median(picl)
        if(zhong<182)
            zhongz=double(zhong)*1.4;
        else
            zhongz=double(zhong);
        end
        %ȡ����ֵ����

        %ָ���任�����Ϊpiclz
        shangbian=zhongz/255;
        xiabian=((zhongz)/255)+((1-zhongz/255)/HUIDUC);
        if(xiabian>1)
            xiabian=1;
        end
        piclz = imadjust(dtemp(:,:,i),[shangbian xiabian],[ ]);
%         figure(13);
%         imshow(piclz);
%         title('ָ���任����ͼ');
        %ָ���任����
       
       
       %��Ŀ�������Ҷ�ֵ�������
       %%
       
       %���õ���䣬ֻ��Ϊ�˵���ʱ��ʾ
       dtemp2= pixel(hang_start:hang_end,lie_start:lie_end,i);
%        figure(3);
%        imshow(dtemp2);
%        title('������Ŀ�������');   
       %���õ���䵽�������

       %�����յĶ�ֵ����ֵ�ٶ�ֵ��һ��
       bw(:, :, i) = im2bw(piclz,ERZHIHUA);              %�Բ�ֵͼ���ֵ�����õ���Ҳ�����������ڴ����ͼ��

       %����
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw(:, :, i),se2);
       %���ͽ���

%        figure(71);
%        imshow(pengzhang);
%        title('���ͺ�');

      %�������������������������򣬽���Щ������1�������ж��������������������ͼ��tf��
      [L,n]=bwlabel(pengzhang);   %����ͨ��nΪ��ͨ��ĸ�����LΪ��bw��С��ͬ�ľ�����ͨ���ѱ��Ϊ��1,2,3....��      
      STATS =regionprops(L,'Area'); %stats��һ������ΪL�Ľṹ���飨ÿλ��һ�����������м�¼������Ӧλ�������С
      aa=struct2cell(STATS);        %�ṹ����ת��ΪԪ�����飨ÿλ��һ������
      bb=cell2mat(aa);              %��Ԫ������ת��Ϊ����һ��n(��ͨ�����)��,��һ�д洢�ľ�����ͨ��1��������ڶ�����ͨ��2����������ơ���
      cc=find((bb>PIXEL_NUM)&(bb<PIXEL_NUM_MAX));         %�ҳ��������Ҫ�����ͨ�򣬷��ص�������Ҫ��������bb�е�λ��,���λ�ö�Ӧ������ʵ������ͨ��ı��
      %��ֵ������

        tf=ismember(L,cc);  %�ж�L�е�Ԫ����û����cc�г��֣�tf��index���Ǻ�Lͬ����С�ľ���
        %tf��Ϊ��ͬλ��1����ͬ��0��������ʵ���ǽ������������ľ�����0��ֻ�������������ľ���
% %         figure(4);
% %         imshow(tf);
% %         title('�ڻ�����Ŀ�����ڶ�ֵ���Ľ��');        
       %%
        %�����洦�����ͼ��tf�����ҳ��������������ͨ����Ϊ��Ŀ���ƶ���ɵģ�
        [L1,n]=bwlabel(tf);
        STATS1 =regionprops(L1,'Area');
        aa1=struct2cell(STATS1);
        bb1=cell2mat(aa1);
        [max_1 max_loction]=max(bb1); 
        [tf1 index]=ismember(L1,max_loction);
         
         
%         figure(5);
%         imshow(L1);
%         title('ȥ�������Сֵ��ʣ������');  
      
    else 
       %%
        d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
        d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
        d(:, :, i)=d(:, :, i)+d2(:, :, i);
        
        
%         figure(2);
%         imshow(d(:, :, i));
%         title('��ֵͼ');

        %��ȷ׷�ٸĽ���   
        %����׷�ٷ�Χ,����d�����bw2
        bw2=d(hangzhong-ZHUIZONG_KUAN:hangzhong+ZHUIZONG_KUAN,liezhong-ZHUIZONG_KUAN :liezhong+ZHUIZONG_KUAN,i);
        %����׷�ٷ�Χ����
       
        
      %%
       %��Ŀ�������Ҷ�ֵ���죬����ͼ��bw2�����ͼ��piclz2
       %��ͼ��ת����һ��
        picl=bw2(:);
%         figure(64);
%         imshow(picl);
%         title('һ��');
        %ͼ��ת����һ�н���

        %ȡ����ֵ
        zhong=median(picl)
        if(zhong<211)
            zhongz=double(zhong)*1.2;
        else
            zhongz=double(zhong);
        end
        %ȡ����ֵ����

        %ָ���任�����Ϊpiclz
        shangbian=zhongz/255;
        xiabian=((zhongz)/255)+((1-zhongz/255)/HUIDUC);
        if(xiabian>1)
            xiabian=1;
        end
        piclz2 = imadjust(bw2,[shangbian xiabian],[]);
%         figure(13);
%         imshow(piclz2);
%         title('ָ���任����ͼ');
       %ָ���任����
       %��Ŀ�������Ҷ�ֵ�������
       %%

        bw3 = im2bw(piclz2,ERZHIHUA); 
         
%         figure(19);
%         imshow(bw3);
%         title('��ֵ�����������');
        
       %����
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw3,se2);
       %���ͽ���

%        figure(71);
%        imshow(pengzhang);
%        title('���ͺ�');
        
        
        %����ͨ��
        [L,n]=bwlabel(pengzhang);   %�����ϴ�׷�٣�Ŀ��λ�õ�ZHUIZONG_KUAN*2�ķ���Ϊ׷����
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>PIXEL_NUM&bb<PIXEL_NUM_MAX);
        tf=ismember(L,cc);
        %����ͨ�򣬻���׷�ٷ�Χ����
        
%         figure(3);
%         imshow(L);
%         title('������Ŀ�������');
%         
%         figure(4);
%         imshow(tf);
%         title('������ͨ��');
        
        [L1,n]=bwlabel(tf);
        STATS1 =regionprops(L1,'Area');
        aa1=struct2cell(STATS1);
        bb1=cell2mat(aa1);
       [max_1 max_loction]=max(bb1);  
    end
    %%
%    �ж��ǲ�������ͨ����������еĻ���������ͨ�����ĳ���ʵ����
    if length(max_loction)   %�ж��ǲ����з�����������ͨ����
         tf1=ismember(L1,max_loction);   %�������ͨ������λ���ҳ���

%%
        %Ѱ������,����Ϊͼ�����Ψһ�����tf1
        [L4,n]=bwlabel(tf1);
        STATS =regionprops(L4,'Centroid');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        liezhongy=round(bb(1,1));
        hangzhongx=round(bb(1,2));
        %Ѱ�����Ľ��������Ϊ���ģ�hangzhongx,liezhongy��
%                 
%         figure(20);
%         imshow(L4);
%         title('�����ͨ����');
%         text(liezhongy,hangzhongx,'+','color','r');
        %��ԭ��ʵ����
        if(i==1)
            hangzhong=hangzhongx+hang_start;
            liezhong=liezhongy+lie_start;
        else
            hangzhong=hangzhongx+hangzhong-ZHUIZONG_KUAN-1;
            liezhong=liezhongy+liezhong-ZHUIZONG_KUAN-1;
        end
        %��ԭ��ʵ�������
        %���Դ��룬ɾ
%         figure(26);
%         imshow(pixels(:, :, :,i));
%         title('ԭͼ');
%         text(liezhong,hangzhong,'+','color','r');
        %���Դ������
    end
    %%  
 
    [l,m]=bwlabel(tf1);
    status=regionprops(l,'BoundingBox');
    
    
    %���Դ��룬ɾ
%     figure;
%     imshow(tf1);
%     title('��ѡ����')
%     rectangle('position',status.BoundingBox,'edgecolor','r');
    %���Դ������
    
    
    
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
    
    
    figure(1);     %��ͼ
    imshow(pixel(:, : ,i));
%     hold on
    rectangle('position',[huitu_lie huitu_hang huitu_kuan huitu_chang],'edgecolor','r');
%     text(huitu_lie+huitu_kuan/2,huitu_hang+huitu_chang/2,'+','color','r');
    
%     figure(1);     %��ͼ
%     imshow(pixels(:, :, :,i));
%     hold on
%     rectangle('Position',[lei_min_last1 hang_min_last1 leikaundu hangkuandu], 'EdgeColor', 'r', 'LineWidth', 2);   %�ӵ�lei_min-5 hang_min-8��ʼ�����Σ���Ϊ20����Ϊ20
% %     rectangle('Position',[lei_min hang_min hang_max-hang_min+20 lei_max-lei_min+20], 'EdgeColor', 'r', 'LineWidth', 2);
    a=i;
    
end


% d(:,:,5)=(abs(pixel(:,:,5)-pixel(:,:,1)));
% bw(:,:,5)=im2bw(d(:,:,5),0.13);
% imshow(bw(:,:,5))