clear all;          %������б���
clc;                %��������
close all           %�ر����д���
avi=VideoReader('W:\KONGYU\se_2\test_file\�·�����_��ʽ_��.mp4');    %ͨ������
%avi=VideoReader('W:\KONGYU\se_2\test_file\����׷��_��.mp4');

%�Ҷ������ϱ߽�����ֵԽ����Ҷ�����ķ�ΧԽխ
HUIDUC=3;
%����
MIN_PIX=120;  %�趨���ڶԽ�������ع��ٵĵ��˳�����ֵ��ע�������˳����Ǳղ��������ͺ��ͼ��
TAG_NUM=20;  %�趨׷��Ŀ�����
BGN=10;      %������ȡ������ͼƬ��
REAL_TAG_NUM=0;   %��ͼ��������շֽ����Ŀ�����������ֱ�����
tag_num =TAG_NUM;  %�м����
ERZHIHUA=0.5;
%%
%�ֽ�ÿһ֡�洢��pixels����ȡ������������֡��
for i=1:avi.NumberOfFrames-100  %ѭ����ȡÿһ֡
    img=read(avi,i);           %��ȡ��ǰ֡
    pixels(:,:,:,i)=img;       %����ǰ֡�洢����ά�����е���Ӧ�㣨��i�㣩ȥ
end

nFrames=size(pixels,4);        %���ؾ���ĵ���ά��size��֡��
rows=size(pixels,1);           %���ؾ���ĵ�һά��size����
cols=size(pixels,2);           %���ؾ���ĵڶ�ά��size����

for i=1:avi.NumberOfFrames-100
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i)));  %RGBת�Ҷ�ͼ�񣬴���pixel
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
        title('����')
        
        figure(7);
        imshow(d2(:,:));
        title('��ȥ����')
    %     d2(:, :) = imsubtract(pixel(:,:,i+1),BACG3(:,:));
    %     figure(8);
    %     imshow(d2(:,:));
    %     title('������2')
    %     d3(:, :) = imsubtract(d2(:, :) , d1(:, :));
    %     figure(9);
    %     imshow(d3(:,:));
    %     title('���')
    %     d4(:, :) = uint8((double(d2(:, :))+double(d1(:, :)))/2);
    %     figure(10);
    %     imshow(d4(:,:));
    %     title('���');
    %     d10(:,:)=pixel(:,:,i)-pixel(:,:,i+1);
    %     figure(17);
    %     imshow(d10(:,:));
    % %     title('���');
    %     d4(:,:)=d4(:, :)+8*d10(:,:);
    %      figure(18);
    %     imshow(d4(:,:));
    %     title('���');
        %%
        %�˲�
        d5=filter2(fspecial('average',4),d2)/255;          %����5*5ģ��ƽ���˲�
        figure(11);
        imshow(d5);
        title('�˲���')
    %%
    %�Ҷ�����
     %��Ŀ�������Ҷ�ֵ���죬����ͼ��d5�����ͼ��piclz
           %��ͼ��ת����һ��
            picl=d5(:);
    %         figure(60);
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
            piclz = imadjust(d5,[shangbian xiabian],[]);
            figure(13);
            imshow(piclz);
            title('ָ���任����ͼ');
            %ָ���任����   
    %�Ҷ��������

     %��ֵ�����Ŷ�
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
       %�жϿ��ŶȽ���
        %��ֵ��
        d6 = im2bw(piclz, ERZHIHUA);
        figure(12);
        imshow(d6);
        title('��ֵ��')
        %%
        %��Ե��ȡ
    %     bianyuantu=edge(d6,'prewitt');%��prewitt���ӽ��б�Ե���
    %     figure;
    %     imshow(bianyuan);
    %     title('��Ե')
    %     d6(1:100,:)=0;
        se=[0,0,1,0,0;0,1,1,1,0;1,1,1,1,1;0,1,1,1,0;0,0,1,0,0];
        pengzhang = imclose(d6,se); %�ղ���
        se2 = strel('disk',3);
        pengzhang2 = imdilate(pengzhang,se2);  %����
    %     
    %     pengzhang(1,:)=0;
    %     pengzhang(rows,:)=0;
    %     pengzhang(:,1)=0;
    %     pengzhang(:,cols)=0;
        figure(13);
        imshow(pengzhang2);
        title('���ͺ�');

       %%
       %���ص���Ŷ�
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
       %�жϿ��ŶȽ���

       %%
        %����ͨ��
         figure(5);
         imshow(pixel(:,:,i));
         [L,n]=bwlabel(pengzhang2);   %����ͨ��nΪ��ͨ��ĸ�����LΪ��bw��С��ͬ�ľ�����ͨ���ѱ��Ϊ��1,2,3....��

    %          figure(22);
    %          imshow(L)
    %          title('������ͨ��')
            if(n~=0)
                STATS =regionprops(L,'Area'); %stats��һ������ΪL�Ľṹ���飨ÿλ��һ�����������м�¼������Ӧλ�������С
                aa=struct2cell(STATS);        %�ṹ����ת��ΪԪ�����飨ÿλ��һ������
                bb=cell2mat(aa);              %��Ԫ������ת��Ϊ����һ��n(��ͨ�����)��,��һ�д洢�ľ�����ͨ��1��������ڶ�����ͨ��2����������ơ���
                cc=find(bb>MIN_PIX);         %�ҳ��������Ҫ�����ͨ�򣬷��ص�������Ҫ��������bb�е�λ��,���λ�ö�Ӧ������ʵ������ͨ��ı��
        %          if length(cc)
                if(cc~=0)
                    tf=ismember(L,cc);  %�ж�L�е�Ԫ����û����cc�г��֣�tf�Ǻ�Lͬ����С�ľ���
            %           %tf��Ϊ��ͬλ��1����ͬ��0��������ʵ���ǽ������������ľ�����0��ֻ�������������ľ���
    %                  figure(21);
    %                  imshow(tf);
    %                  title('ȥ�������С�������ͨ��');
                     %�������ж��Ƿ��Ƿ����������������

                     %%
                     %ȡ������Ԥ�Ƶ�n��tag_num����Ŀ������
                    [L1,n]=bwlabel(tf);
                    if(n~=0)
                        STATS1 =regionprops(L1,'Area'); 
                        aa1=struct2cell(STATS1);
                        bb1=cell2mat(aa1);
                        num_can_do=length(bb1);%��Ŀ���޸���1
                        if(num_can_do<tag_num)
                            tag_num=num_can_do;
                        end
                        for tag_num_i=1:tag_num     
                            [max_1(tag_num_i) max_loction]=max(bb1);%max_1�����ֵ��max_location�����ֵ��һ�γ��ֵ�λ��,��ʵ������������ͨ��ı��
                            tf1(:,:,tag_num_i)=ismember(L1,max_loction);%�����������ͨ����1������λ����0����ͼ��tf1��tf1�о�ֻ�У�����Ŀ�������
                            bb1(max_loction)=[1];

                            figure(20);
                            imshow(tf1(:,:,tag_num_i));
                            title('�ֽ����Ŀ��');
                        end
                        tf2=tf1(:,:,1);
                        for tag_num_i=2:tag_num
                            tf2=tf2+tf1(:,:,tag_num_i);
                        end
                        tag_num=TAG_NUM;

                        figure;
                        imshow(tf2);
                        title('����Ŀ��')
                       %%

                         %����Ŀ���
                         [l,m]=bwlabel(tf2,8);
                         status=regionprops(l,'BoundingBox');
                         tar_position_temp=struct2cell(status);
                         tar_position=cell2mat(tar_position_temp);    %(�У��У��п��п�)

                         REAL_TAG_NUM=(length(tar_position))/4;

                         figure(71);
                         imshow(pixel(:,:,i));
                         title('�����');
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

