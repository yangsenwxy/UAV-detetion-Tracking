clear all;          %������б���
clc;                %��������
close all           %�ر����д���
avi=VideoReader('D:\KONGYU\se_2\yiceshi\cs0409.mp4');    %��ȡ��Ƶ��avi
t1=zeros(350,1); 
t4=t1;
KUAN_4_CHA1=zeros(350,4);
MAXKUANG_4MJC1=zeros(350,2);
KUAN_4_MJC1=zeros(350,2);
TAR_LOS_TIMES=0;
%�ֶ�����ı߽�,��װ�������������
% hang_start=274;   %�·�����_��ʽ_��,ȥ��ǰ50֡��ͨ��
% hang_end=301;
% lie_start=248;                                                                                                                                                                                                                                                                                    00;
% % lie_end=283;
% 
% hang_start=204;%������_�̶���_��ʽ��_1��50��ȥ��ǰ50֡,ͨ��,
% hang_end=228;
% lie_start=307;                                                                                                                                                                                                                                                                                    00;
% lie_end=336;

% hang_start=308;%1������_��ʽ��20��ȥ��ǰ50֡��ͨ��
% hang_end=339;
% lie_start=337;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;
%����


hang_start=350;  %20170329_��������ȥ��ǰ50֡��ͨ��
hang_end=390;
lie_start=225;                                                                                                                                                                                                                                                                                    00;
lie_end=265;

% hang_start=236;  %20170329_������2��ȥ��ǰ50֡��ͨ��
% hang_end=253;
% lie_start=448;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;

% 
% hang_start=341;  %20170329_������2cuo��ȥ��ǰ50֡��ͨ��
% hang_end=383;
% lie_start=414;                                                                                                                                                                                                                                                                                    00;
% lie_end=445;

% hang_start=349;  %20170401_clip��ȥ��ǰ50֡��ͨ��
% hang_end=375;
% lie_start=420;                                                                                                                                                                                                                                                                                    00;
% lie_end=442;

%�Ҷ������ϱ߽�����ֵԽ����Ҷ�����ķ�ΧԽխ
HUIDUC=2;
%����


%׷�����Ŀ�ȣ������ֶ����ƿ�Ŀ����ȷ��
KUANG_2=round(max((hang_end-hang_start),(lie_end-lie_start))*0.5);

% %��ȷ׷�ٿ�Ŀ��
JINGQUE_KUAN=round(KUANG_2*0.5);

K_JINGQUE=KUANG_2-JINGQUE_KUAN;

%���Ͳ���
% PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/20);
PENGZHANG=3;
%���ͽ���

PIXEL_NUM=60;  %ʶ��Ŀ����������ص����
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*0.5);%ʶ��Ŀ���������ص����



ERZHIHUA=0.08;  %Ԥ���ֵ������ֵ�����Ǿ���ֵ����1kmĿ��ʶ��������Ч�ģ����ظ���
ERZHIHUA2=0.8;

YOUXIAO=1;   %Ŀ����Ч
TAR_LOST=0;  %Ŀ��δ��ʧ
TAR_LOST_TIMES=0;   %Ŀ�궪ʧ֡��ͳ��
TAR_LOST_SIGN=5;  %ȷ��Ŀ���Ѷ�ʧ��֡�����ޣ�����TAR_LOST_SIGN֡Ŀ�����ʧʱ����ΪĿ��ȷʵ�Ѿ���ʧ�������ⲿ����

TAR_ON_HOLD=0;

QUEREN=0;

%��ȷ׷������ȱ���������Ӧ����
%%
%�ֽ�ÿһ֡�洢��pixels����ȡ������������֡��
for i=1:avi.NumberOfFrames-50  %ѭ����ȡÿһ֡
    img=read(avi,i);           %��ȡ��ǰ֡
    pixels(:,:,:,i)=img;       %����ǰ֡�洢����ά�����е���Ӧ�㣨��i�㣩ȥ
end

nFrames=size(pixels,4);        %���ؾ���ĵ���ά��size��֡��
rows=size(pixels,1);           %���ؾ���ĵ�һά��size����
cols=size(pixels,2);           %���ؾ���ĵڶ�ά��size����

for i=1:nFrames-50
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i+50)));  %RGBת�Ҷ�ͼ�񣬴���pixel
end

k=1;

%%
%��ѭ����������������Ƶ������ѭ��
for i =1 : nFrames-50     
   %%
    if i==1      %�ж��Ƿ�Ϊ��һ֡ͼ��
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);

%          figure(2);
%          imshow(d(:, :, i));
%          title('��ֵͼ1');
%          
       %%
       %����Ŀ�������
       dtemp(:,:,i)=d(hang_start:hang_end,lie_start:lie_end,i);
%        figure(3);
%        imshow(dtemp(:, :, i));
%        title('Ŀ������ֵͼ');
       %����Ŀ��������
       %%
       %��Ŀ�������Ҷ�ֵ���죬����ͼ��dtemp(:,:,i)�����ͼ��piclz
       %��ͼ��ת����һ��
        picl=dtemp(:);
        SUB_GRAY_AVR_LAST=sum(picl)/((hang_end-hang_start)*(lie_end-lie_start));
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
        piclz = imadjust(dtemp(:,:,i),[shangbian xiabian],[]);
        
        
        
        %ָ���任����
       %��Ŀ�������Ҷ�ֵ�������
       %%
       %�����յĶ�ֵ����ֵ�ٶ�ֵ��һ��
       bw(:, :, i) = im2bw(piclz,ERZHIHUA);              %�Բ�ֵͼ���ֵ�����õ���Ҳ�����������ڴ����ͼ��

       %����
       se2 = strel('disk',PENGZHANG);
       pengzhang = imdilate(bw(:, :, i),se2);
       %���ͽ���

      %�������������������������򣬽���Щ������1�������ж��������������������ͼ��tf��
      [L,n]=bwlabel(pengzhang);   %����ͨ��nΪ��ͨ��ĸ�����LΪ��bw��С��ͬ�ľ�����ͨ���ѱ��Ϊ��1,2,3....��      
      STATS =regionprops(L,'Area'); %stats��һ������ΪL�Ľṹ���飨ÿλ��һ�����������м�¼������Ӧλ�������С
      aa=struct2cell(STATS);        %�ṹ����ת��ΪԪ�����飨ÿλ��һ������
      bb=cell2mat(aa);              %��Ԫ������ת��Ϊ����һ��n(��ͨ�����)��,��һ�д洢�ľ�����ͨ��1��������ڶ�����ͨ��2����������ơ���
      cc=find((bb>PIXEL_NUM)&(bb<PIXEL_NUM_MAX));         %�ҳ��������Ҫ�����ͨ�򣬷��ص�������Ҫ��������bb�е�λ��,���λ�ö�Ӧ������ʵ������ͨ��ı��
      %��ֵ������

        tf=ismember(L,cc);  %�ж�L�е�Ԫ����û����cc�г��֣�tf��index���Ǻ�Lͬ����С�ľ���
        %tf��Ϊ��ͬλ��1����ͬ��0��������ʵ���ǽ������������ľ�����0��ֻ�������������ľ���  
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
        YOUXIAO=1;
      
    else 
       %%
            d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
            d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
            d(:, :, i)=d(:, :, i)+d2(:, :, i);

%             figure(5);
%             imshow(d(:, :, i));
%             title('��ֵͼ');



           %�������Ȩ��  

            %��������Ȩ����,����d�����bw4
            bw4=d(hangzhong-JINGQUE_KUAN:hangzhong+JINGQUE_KUAN,liezhong-JINGQUE_KUAN :liezhong+JINGQUE_KUAN,i);
            %����׷�ٷ�Χ����

            
            figure(6);
            imshow(bw4);
            title('����Ȩ����');


           %��Ŀ�������Ҷ�ֵ���죬����ͼ��bw2�����ͼ��piclz2
           %��ͼ��ת����һ��
            picl4=bw4(:);
            %ȡ����ֵ
            zhong4=median(picl4)
            if(zhong4<211)
                zhongz4=double(zhong4)*1.2;
            else
                zhongz4=double(zhong4);
            end
            %ȡ����ֵ����

            %ָ���任�����Ϊpiclz
            shangbian4=zhongz4/255;
            xiabian4=((zhongz4)/255)+((1-zhongz4/255)/HUIDUC);
            if(xiabian4>1)
                xiabian4=1;
            end
            piclz4 = imadjust(bw4,[shangbian4 xiabian4],[]);
            
            
            figure(7);
            imshow(piclz4);
            title('����Ȩ�����Ҷ�����');
            
            
            picl14=piclz4(:);
            %ȡ����ֵ
            zhong14=median(picl14)
            if(zhong14<211)
                zhongz14=double(zhong14)*1.2;
            else
                zhongz14=double(zhong14);
            end
            %ȡ����ֵ����

            %ָ���任�����Ϊpiclz
            shangbian14=zhongz14/255;
            xiabian14=((zhongz14)/255)+((1-zhongz14/255)/HUIDUC);
            if(xiabian14>1)
                xiabian14=1;
            end
            piclz14 = imadjust(piclz4,[shangbian4 xiabian14],[]);

            %�������Ȩ�޽���


            figure(7);
            imshow(piclz14);
            title('����Ȩ�������λҶ�����');

            %����׷�ٷ�Χ,����d�����bw2
            bw2=d(hangzhong-KUANG_2:hangzhong+KUANG_2,liezhong-KUANG_2 :liezhong+KUANG_2,i);
            bw5=bw2;
            
    %         bw7(K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1),K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN+1))=bw7+piclz4;
            for z=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                for y=K_JINGQUE:(K_JINGQUE+2*JINGQUE_KUAN)
                    bw5(z,y)=piclz14(z-K_JINGQUE+1,y-K_JINGQUE+1);
                end
            end
            bw6=bw5;
            %����׷�ٷ�Χ����


            figure(8);
            imshow(bw6);
            title('׷��������Ȩ��������');
            
            %s1,����ͼ4�ſ���ҶȾ�ֵgray_sum;
            tongjiqu=pixel(HUITU_HANG:(HUITU_HANG+HUITU_CHANG),HUITU_LIE:(HUITU_LIE+HUITU_KUAN),i);
            gray_avr=sum(tongjiqu(:))/(HUITU_KUAN*HUITU_CHANG);


          %%
           %��Ŀ�������Ҷ�ֵ���죬����ͼ��bw2�����ͼ��piclz2
           %��ͼ��ת����һ��
            picl=bw2(:);
            picl6=bw6(:);
            %ͼ��ת����һ�н���

            %s2,��׷�ٿ�����ҶȾ�ֵsub_gray_sum
            sub_gray_avr=sum(picl)/((KUANG_2*2+1)*(KUANG_2*2+1));
            %s3,�ж�׷�����ҶȾ�ֵ���ϴεĻҶȾ�ֵ������
%xw             
                t=abs(SUB_GRAY_AVR_LAST-sub_gray_avr);
                t1(i,:)=t;
                t2=median(t1);
                t3=max(t1);
                t4(i,:)=t3;
                t5=sum(t4(:))/i;
             
            if ((t>t2)&&(t<t5))%׷���������ҶȾ�ֵ���
%                if (t<t5)%׷���������ҶȾ�ֵ���
                %�������˶�
                SUB_GRAY_AVR_LAST=sub_gray_avr;
                YOUXIAO=1;%�Ƿ�Ŀ������ֵͼ���������������㣨����ΪĿ��ȷʵ�������㹻���ƶ�������ִ�к�������
                TAR_ON_HOLD=0;
                %ȡ����ֵ
                zhong=median(picl6);
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
                piclz2 = imadjust(bw6,[shangbian xiabian],[]);
               %ָ���任����
               %��Ŀ�������Ҷ�ֵ�������
%%

                bw3 = im2bw(piclz2,ERZHIHUA); 
                
                
                figure(9);
                imshow(bw3);
                title('����Ȩ�������λҶ�������ֵ��');
                
                
                %����
                se2 = strel('disk',PENGZHANG);
                pengzhang = imdilate(bw3,se2);
                %���ͽ���
                
                figure(10);
                imshow(pengzhang);
                title('����Ȩ�������λҶ�������ֵ��������');

                %����ͨ��
                [L,n]=bwlabel(pengzhang);   %�����ϴ�׷�٣�Ŀ��λ�õ�KUANG_2*2�ķ���Ϊ׷����
                STATS =regionprops(L,'Area');
                aa=struct2cell(STATS);
                bb=cell2mat(aa);
                cc=find(bb>PIXEL_NUM&bb<PIXEL_NUM_MAX);
                tf=ismember(L,cc);
                %����ͨ�򣬻���׷�ٷ�Χ����
                 if length(cc)==0
                     YOUXIAO=0;
                 else
                    [L1,n]=bwlabel(tf);
                    STATS1 =regionprops(L1,'Area');
                    aa1=struct2cell(STATS1);
                    bb1=cell2mat(aa1);
                    [max_1 max_loction]=max(bb1);  
                 end

            else%��ǰ���ܾ�ֹ��ʧ
                YOUXIAO=0;
                TAR_ON_HOLD=1;
            end
%             if YOUXIAO==0
%                
%             end
    end
    %%
%    �ж��ǲ�������ͨ����������еĻ���������ͨ�����ĳ���ʵ����
    if (YOUXIAO==1)&&(length(max_loction))   %�ж��ǲ����з�����������ͨ����
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
        
        %��ԭ��ʵ����
        if(i==1)
            hangzhong=hangzhongx+hang_start;
            liezhong=liezhongy+lie_start;
             
        else
            hangzhong=hangzhongx+hangzhong-KUANG_2-1;
            liezhong=liezhongy+liezhong-KUANG_2-1;
        end
        %��ԭ��ʵ�������

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
               TAR_LOST=0;    %Ŀ���Ƿ�ʧ����ʧ��1
               GRAY_AVR=gray_avr;
               TAR_ON_HOLD=1;   %������������һ�η���Ŀ��onhold������һ����ʼ״̬Ϊ0����ʾû�о�ֹ����Ϊ2������������ٵ�2�η���Ŀ��onhold
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
    figure(1);     %��ͼ
    imshow(pixel(:, : ,i));
    fprintf('i= %d\n',i);
    rectangle('position',[HUITU_LIE,HUITU_HANG HUITU_KUAN HUITU_CHANG],'edgecolor','r');
    
    
    
    %���ԣ�ɾ
%     figure(2);     %��ͼ
%     imshow(d(:, : ,i));
%     rectangle('position',[HUITU_LIE HUITU_HANG HUITU_KUAN HUITU_CHANG],'edgecolor','r');
    %���Խ���
    
 
%     YOUXIAO=0;
    if i>58
       a=1;                                        
   end
end