clear all;          %������б���
clc;                %��������
close all           %�ر����д���
avi=mmreader('W:\KONGYU\se_2\test_file\20170401_clip.mp4');    %��ȡ��Ƶ��avi

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
% lie_start=337;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;
%����


% hang_start=323;  %20170329_��������ȥ��ǰ50֡��ͨ��
% hang_end=342;
% lie_start=363;                                                                                                                                                                                                                                                                                    00;
% lie_end=395;

% hang_start=236;  %20170329_������2��ȥ��ǰ50֡��ͨ��
% hang_end=253;
% lie_start=448;                                                                                                                                                                                                                                                                                    00;
% lie_end=465;


hang_start=341;  %20170329_������2��ȥ��ǰ50֡��ͨ��
hang_end=383;
lie_start=414;                                                                                                                                                                                                                                                                                    00;
lie_end=445;

%�Ҷ������ϱ߽�����ֵԽ����Ҷ�����ķ�ΧԽխ
HUIDUC=2;
%����


%׷�����Ŀ�ȣ������ֶ����ƿ�Ŀ����ȷ��
ZHUIZONG_KUAN=round(max((hang_end-hang_start),(lie_end-lie_start))*0.75);

% %��ȷ׷�ٿ�Ŀ��
% JINGQUE_KUAN=round(ZHUIZONG_KUAN*1.25);

%���Ͳ���
PENGZHANG=round((max((hang_end-hang_start),(lie_end-lie_start)))/10);
%���ͽ���

PIXEL_NUM=2;  %ʶ��Ŀ����������ص����
PIXEL_NUM_MAX=round(((((hang_end-hang_start)*(lie_end-lie_start))))*0.9);%ʶ��Ŀ���������ص����


ERZHIHUA=0.1;  %Ԥ���ֵ������ֵ�����Ǿ���ֵ����1kmĿ��ʶ��������Ч�ģ����ظ���

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



for i =1 : nFrames-100      
         d(:, :, i) = (abs(pixel(:,:,i+1) - pixel(:,:,i)));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
         d2(:, :, i) = (abs(pixel(:,:,i) - pixel(:,:,i+1)));
         d(:, :, i)=d(:, :, i)+d2(:, :, i);
%        d100(:, :, i) = imsubtract(pixel(:,:,i+1),pixel(:,:,i));   %��ǰ֡����һ֡�ĲdΪ����Ĳ�ֵͼ��
%        d1(:,:,i)=(abs(pixel(:,:,i) - pixel(:,:,i+1)));
%        d3(:,:,i)=uint8(((double(d(:, :, i))+double(d1(:, :, i))))/2);

         figure(2);
         imshow(d(:, :, i));
         title('��ֵͼ1');
end