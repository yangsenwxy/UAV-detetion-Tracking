clear all;          %������б���
clc;                %��������
close all           %�ر����д���
avi=VideoReader('W:\KONGYU\se_2\test_file\������_�̶���_��ʽ��_1.mp4');    %��ȡ��Ƶ��avi

ERZHIHUA=0.05;
JINGQUE_L=120;
JINGQUE_2L=JINGQUE_L*2+1;
%%
%�ֽ�ÿһ֡�洢��pixels����ȡ������������֡��
for i=1:avi.NumberOfFrames-200  %ѭ����ȡÿһ֡
    img=read(avi,i);           %��ȡ��ǰ֡
    pixels(:,:,:,i)=img;       %����ǰ֡�洢����ά�����е���Ӧ�㣨��i�㣩ȥ
end

nFrames=size(pixels,4);        %���ؾ���ĵ���ά��size��֡��
rows=size(pixels,1);           %���ؾ���ĵ�һά��size����
cols=size(pixels,2);           %���ؾ���ĵڶ�ά��size����

for i=1:nFrames
  pixel(:,:,i)=(rgb2gray(pixels(:,:,:,i)));  %RGBת�Ҷ�ͼ�񣬴���pixel
%   figure;
%   imshow(pixels(:,:,:,1));
%   title('yuantu');
%   figure;
%   imshow(pixel(:,:,i));
%   title('huidu')
end

k=1;

%%
%��ѭ����������������Ƶ������ѭ��
for i = 1 : nFrames-4       
    d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
    d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
    d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);
    bw(:, :, i) = im2bw(d(:, :, i), ERZHIHUA);                %�Բ�ֵͼ���ֵ�����õ���Ҳ�����������ڴ����ͼ��
%     
%      figure;
%      imshow(d(:,:,i));
%      title('�ҶȲ�ֵͼ1')
     
%      figure;
%      imshow(d1(:,:,i));
%      title('�ҶȲ�ֵͼ2')
%      
%      figure;
%      imshow(d3(:,:,i));
%      title('��ֵͼ���')
%      
% %      imhist(d(:,:,i));
%      figure;
%      imshow(pixel(:,:,i))
%      title('�Ҷ�ͼ1');
%      figure;
%      imshow(pixel(:,:,i+1))
%      title('�Ҷ�ͼ2')
%      figure;
%      imshow(bw(:, :, i));
%      title('��ֵ���Ĳ�ֵͼ');
% %     
   %%
    if i==1      %�ж��Ƿ�Ϊ��һ֡ͼ��
        
        
        %�������������������������򣬽���Щ������1�������ж��������������������ͼ��tf��
        [L,n]=bwlabel(bw(:, :, i));   %����ͨ��nΪ��ͨ��ĸ�����LΪ��bw��С��ͬ�ľ�����ͨ���ѱ��Ϊ��1,2,3....��
        
%          figure;
%          imshow(L)
%          title('������ͨ��')
%         
        STATS =regionprops(L,'Area'); %stats��һ������ΪL�Ľṹ���飨ÿλ��һ�����������м�¼������Ӧλ�������С
        aa=struct2cell(STATS);        %�ṹ����ת��ΪԪ�����飨ÿλ��һ������
        bb=cell2mat(aa);              %��Ԫ������ת��Ϊ����һ��n(��ͨ�����)��,��һ�д洢�ľ�����ͨ��1��������ڶ�����ͨ��2����������ơ���
        cc=find(bb>2&bb<100);         %�ҳ��������Ҫ�����ͨ�򣬷��ص�������Ҫ��������bb�е�λ��,���λ�ö�Ӧ������ʵ������ͨ��ı��
%          if length(cc)
         tf=ismember(L,cc);  %�ж�L�е�Ԫ����û����cc�г��֣�tf��index���Ǻ�Lͬ����С�ľ���
          %tf��Ϊ��ͬλ��1����ͬ��0��������ʵ���ǽ������������ľ�����0��ֻ�������������ľ���
%          figure;
%          imshow(tf);
%          title('ȥ�������С����Ķ�ֵͼ');
         %�������ж��Ƿ��Ƿ����������������
         
         
       %%
        %�����洦�����ͼ��tf�����ҳ��������������ͨ����Ϊ��Ŀ���ƶ���ɵģ�
        [L,n]=bwlabel(bw(:, :, i));
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%          if length(cc)
          [tf index]=ismember(L,cc);
%          figure,
%          imshow(tf);
%          title('����Ԥ��Ŀ��')
         
        se=[0 0 1 0 0,0 0 1 0 0,1 1 1 1 1,0 0 1 0 0,0 0 1 0 0];    %�ṹ��Ԫ
        pengzhang = imclose(tf,se);     %��se
        se2 = strel('disk',10);
        tf= imdilate(tf,se2);
%         figure;
%         imshow(tf);
%         title('�ղ���');
         
         
        [L1,n]=bwlabel(tf);
         STATS1 =regionprops(L1,'Area');
         aa1=struct2cell(STATS1);
         bb1=cell2mat(aa1);
         [max_1 max_loction]=max(bb1); 
         [tf1 index]=ismember(L1,max_loction);
        
      
    else 
       %%
        %����ͨ�򣬱궨���¿��ܵ�����ͬ��
        [L,n]=bwlabel(bw(:, :, i));
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%       if length(cc)
        tf=ismember(L,cc);
       
        if(((hang_min-20)<1)||((lei_min-20)<1)||((hang_min+20)>rows)||((lei_min+20)>cols))
            msgbox('�ɳ���Χ');
        end
        [L,n]=bwlabel(tf(hang_min-20:hang_min+20,lei_min-20 :lei_min+20));   %hang_min��lei_min���ϴλ��Ĳο��㣬����ʵ���꣬��ʱ�����ϴβο���Ϊ����41*41�ķ�Χ��Ѱ��Ŀ��
        STATS =regionprops(L,'Area');
        aa=struct2cell(STATS);
        bb=cell2mat(aa);
        cc=find(bb>2&bb<100);
%       if length(cc)
        tf=ismember(L,cc);
%         figure,
%         imshow(tf);
%         title('����Ԥ��Ŀ��');
        
       se=[0 0 1 0 0,0 0 1 0 0,1 1 1 1 1,0 0 1 0 0,0 0 1 0 0];    %�ṹ��Ԫ
        pengzhang = imclose(tf,se);     %��se
        se2 = strel('disk',10);
        tf= imdilate(tf,se2);
%         figure;
%         imshow(tf);
%         title('�ղ���');
        
%       ȷ���Ƿ���ڷ���Ԥ�ڵ�Ŀ��
       [L1,n]=bwlabel(tf);
       STATS1 =regionprops(L1,'Area');
       aa1=struct2cell(STATS1);
       bb1=cell2mat(aa1);
      [max_1 max_loction]=max(bb1);  
    end
    %%
    %�ж��ǲ�������ͨ����������еĻ���������ͨ�����ĳ���ʵ����
    if length(max_loction)   %�ж��ǲ����з�����������ͨ����
         tf1=ismember(L1,max_loction);   %�������ͨ������λ���ҳ���

%%
          %ȡ�����ϵĵ�Ϊ�ο���
          if i==1
              gg=find(tf1==1)';   %�ҳ�����������λ�ñ�ţ������һά����gg�У�����gg���к��м�������С���ܴﵽcow*lei
              lei=fix(gg./rows)+1; %��0ȡ��
              lei_min=min(lei);   %ȡ����С��ֵ
%               lei_min_1=lei_min;
              lei_max=max(lei);
              hang=rem(gg,rows)+1; %ȡ�࣬��С����ֵ
              hang_min=min(hang);
%               hang_min_1=hang_min;
              hang_max=max(hang);
          else 
              gg=find(tf1==1)';
              lei=fix(gg./41)+1-21+lei_min_last; %����ʵ�ο�����
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
    %��ȷ׷�ٸĽ���
    %mubiaoqu=d3(hang_min_last-100:hang_min_last+100,lei_min_last-100:lei_min_last+100,i);
     if(((hang_min_last-JINGQUE_L)<1)||((lei_min_last-JINGQUE_L)<1)||((hang_min_last+JINGQUE_L)>rows)||((lei_min_last+JINGQUE_L)>cols))
            msgbox('�ɳ���Χ');
     end
      edge_top=hang_min_last-JINGQUE_L;
      edge_bottom=hang_min_last+JINGQUE_L;
      edge_left=lei_min_last-JINGQUE_L;
      edge_right=lei_min_last+JINGQUE_L;
    %%
    %��ʼ����Ŀ��
    mubiaoqu=bw(edge_top:edge_bottom,edge_left:edge_right,i);
%     figure;
%     imshow(mubiaoqu);
%     title('Ŀ����')
    bianyuantu=edge(mubiaoqu,'prewitt');%��prewitt���ӽ��б�Ե���
    
%     figure(3);
%     imshow(bianyuantu);
%     title('��Ե')
%     
    se=[0 0 0 0 1 0 0 0 0,0 0 0 0 1 0 0 0 0,0 0 0 1 1 1 0 0 0,0 0 1 1 1 1 1 0 0,1 1 1 1 1 1 1 1 1,0 0 1 1 1 1 1 0 0,0 0 0 1 1 1 0 0 0,0 0 0 0 1 0 0 0 0,0 0 0 0 1 0 0 0 0];    %�ṹ��Ԫ
    bicaozuo = imclose(bianyuantu,se);     %��se
    
%     figure;
%     imshow(bicaozuo);
%     title('�ղ���')
    
    se2 = strel('disk',20);
    pengzhang = imdilate(bicaozuo,se2);  %����
%     figure(2);
%     imshow(pengzhang);
%     title('���ͺ�');
%     
% 
    pengzhang(1,:)=0;
    pengzhang(JINGQUE_2L,:)=0;
    pengzhang(:,1)=0;
    pengzhang(:,JINGQUE_2L)=0;

    %���������
    [PL1,n]=bwlabel(pengzhang);
    STATS2 =regionprops(PL1,'Area'); 
    paa1=struct2cell(STATS2);
    pbb1=cell2mat(paa1);
    cc=find(pbb1>200&pbb1<10000);
    tf3=ismember(PL1,cc);
%     figure(3);
%     imshow(tf3);
%     title('�����');
%     %��������
    [L4,n]=bwlabel(tf3);
    STATS4 =regionprops(L4,'Area');
    aa4=struct2cell(STATS4);
    bb4=cell2mat(aa4);
    [max_val max_loc1]=max(bb4);
    tf4=ismember(L4,max_loc1);
%     figure(4);
%     imshow(tf4);
%     title('�����');
    %���ʵ������ֵ
    weiz=find(tf4==1)';   %�ҳ�����������λ�ñ�ţ������һά����gg�У�����gg���к��м�������С���ܴﵽcow*lei
    lei=fix(weiz./JINGQUE_2L)+1; %��0ȡ��
    lei_min1=min(lei)-JINGQUE_L+lei_min;   %ȡ����С��ֵ
    lei_max1=max(lei)+lei_max-JINGQUE_L;
    hang=rem(weiz,JINGQUE_2L)+1; %ȡ�࣬��С����ֵ
    hang_min1=min(hang)+hang_min-JINGQUE_L;
    hang_max1=max(hang)+hang_max-JINGQUE_L;
    leikuandu=max(lei)-min(lei);
    hangkuandu=max(hang)-min(hang);
   %% 
    %ȷ���˴��Ƿ��ҵ�Ŀ�꣬���û���ҵ��������ϴ�λ��
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
    figure(1);     %��ͼ
    imshow(pixels(:, :, :, i));
    hold on
    %%
    %�궨ͼ��
    rectangle('Position',[lei_min1 hang_min1 leikuandu hangkuandu], 'EdgeColor', 'r', 'LineWidth', 2);   %�ӵ�lei_min-5 hang_min-8��ʼ�����Σ���Ϊ20����Ϊ20
%     rectangle('Position',[lei_min hang_min hang_max-hang_min+20 lei_max-lei_min+20], 'EdgeColor', 'r', 'LineWidth', 2);
    a1212a=1;
    
end


% d(:,:,5)=(abs(pixel(:,:,5)-pixel(:,:,1)));
% bw(:,:,5)=im2bw(d(:,:,5),0.13);
% imshow(bw(:,:,5))