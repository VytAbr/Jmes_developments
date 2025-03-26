clc; 
clear all; 
close all;
zodziu_sk = 100; s=zeros(128,256); eile = 12;
A2 = zeros(zodziu_sk,128,eile);  C2 = zeros(zodziu_sk,128,eile); 
A2s = zeros(1,zodziu_sk);
A4_1 = zeros(zodziu_sk,128,eile);  
A4s = zeros(1,zodziu_sk);
A4_2 = zeros(zodziu_sk,128,eile);  
A4_3 = zeros(zodziu_sk,128,eile);  
A4_4 = zeros(zodziu_sk,128,eile);  
A4_5 = zeros(zodziu_sk,128,eile);  
A4_6 = zeros(zodziu_sk,128,eile);  
A4_7 = zeros(zodziu_sk,128,eile);  
A4_8 = zeros(zodziu_sk,128,eile);  
A4_9 = zeros(zodziu_sk,128,eile);  
A4_10 = zeros(zodziu_sk,128,eile);  
for j=1:zodziu_sk % 1:100
%     [data_y4 length_y] = pavad_normalization_0215(j,'/M030',3);
    [data_y4_1 length_y] = pavad_normalization_0125(j,'/V020',1);
    [data_y4_2 length_y] = pavad_normalization_0125(j,'/V020',2);
    [data_y4_3 length_y] = pavad_normalization_0125(j,'/V020',3);
    [data_y4_4 length_y] = pavad_normalization_0125(j,'/V020',4);
    [data_y4_5 length_y] = pavad_normalization_0125(j,'/V020',5);
    [data_y4_6 length_y] = pavad_normalization_0125(j,'/V020',6);
    [data_y4_7 length_y] = pavad_normalization_0125(j,'/V020',7);
    [data_y4_8 length_y] = pavad_normalization_0125(j,'/V020',8);
    [data_y4_9 length_y] = pavad_normalization_0125(j,'/V020',9);
        [data_y4_10 length_y] = pavad_normalization_0125(j,'/V02',10);
%     [data_y2 length_y] = pavad_normalization_0125(j,'/V030',7);
    len_y = floor(length_y);
    for i = 1:128 % 2:2 % 1:128
% % % % % % % % -----------------------------------------
        s = fix((data_y4_1((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_1(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_2((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_2(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_3((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_3(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_4((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_4(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_5((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_5(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_6((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_6(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_7((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_7(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_8((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_8(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_9((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_9(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        s = fix((data_y4_10((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        A4_10(j,i,:) = mfcc(s,128); %a4(:,:)=A4_1(j,i,:);
        A4s(j) = len_y;
% % % % % % % % -----------------------------------------
% %         t = fix((data_y2((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
% %         t_tr = 0;
% %         A2(j,i,:) = mfcc(t,128);       a2(:,:)=A2(j,i,:);
% %         A2s(j) = len_y;
    end
% %     disp(j);% commented by AS 2013 03 27
end
% % % figure(7);
% % % plot(t); grid on;

error = zeros(zodziu_sk,zodziu_sk); 
word_nr=zeros(1,zodziu_sk); 
word_nr2=zeros(1,zodziu_sk); 
cnt = 0; cnt2 = 0; rate = 0; rate2 = 0; k = 1; k2 = 1;
cnt_som = 0; rate_som = 0;
% int_error = zeros(zodziu_sk,zodziu_sk); 
% serieA = zeros(128,12); serieB = zeros(128,12);
%% SOM klasifikavimas
klases = [6 1];
% % net = selforgmap([klases 1], 400, klases-1);
net = selforgmap(klases,200,5);
net = train(net,reshape(A4_1,12800,[])');
% % Results_A2 = sim(net, reshape(A2,12800,[])');
% % Results_A4 = sim(net, reshape(A4,12800,[])');
% % [yy,xx]=find(Results_A2);
close all
SOM_A4_1=zeros(128,100);
SOM_A4_2=zeros(128,100);
SOM_A4_3=zeros(128,100);
SOM_A4_4=zeros(128,100);
SOM_A4_5=zeros(128,100);
SOM_A4_6=zeros(128,100);
SOM_A4_7=zeros(128,100);
SOM_A4_8=zeros(128,100);
SOM_A4_9=zeros(128,100);
SOM_A4_10=zeros(128,100);
for zodzio_nr = 1:zodziu_sk
    A4_1_tmp = A4_1(zodzio_nr,:,:);
    Results_A4_1_tmp = sim(net, reshape(A4_1_tmp,128,[])');
    [yy,xx]=find(Results_A4_1_tmp);
%     figure(zodzio_nr),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_1(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_2_tmp = A4_2(zodzio_nr,:,:);
    Results_A4_2_tmp = sim(net, reshape(A4_2_tmp,128,[])');
    [yy,xx]=find(Results_A4_2_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_2(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_3_tmp = A4_3(zodzio_nr,:,:);
    Results_A4_3_tmp = sim(net, reshape(A4_3_tmp,128,[])');
    [yy,xx]=find(Results_A4_3_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_3(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_4_tmp = A4_4(zodzio_nr,:,:);
    Results_A4_4_tmp = sim(net, reshape(A4_4_tmp,128,[])');
    [yy,xx]=find(Results_A4_4_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_4(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_5_tmp = A4_5(zodzio_nr,:,:);
    Results_A4_5_tmp = sim(net, reshape(A4_5_tmp,128,[])');
    [yy,xx]=find(Results_A4_5_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_5(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_6_tmp = A4_6(zodzio_nr,:,:);
    Results_A4_6_tmp = sim(net, reshape(A4_6_tmp,128,[])');
    [yy,xx]=find(Results_A4_6_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_6(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_7_tmp = A4_7(zodzio_nr,:,:);
    Results_A4_7_tmp = sim(net, reshape(A4_7_tmp,128,[])');
    [yy,xx]=find(Results_A4_7_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_7(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_8_tmp = A4_8(zodzio_nr,:,:);
    Results_A4_8_tmp = sim(net, reshape(A4_8_tmp,128,[])');
    [yy,xx]=find(Results_A4_8_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_8(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_9_tmp = A4_9(zodzio_nr,:,:);
    Results_A4_9_tmp = sim(net, reshape(A4_9_tmp,128,[])');
    [yy,xx]=find(Results_A4_9_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_9(:,zodzio_nr)=yy;
end
for zodzio_nr = 1:zodziu_sk
    A4_10_tmp = A4_10(zodzio_nr,:,:);
    Results_A4_10_tmp = sim(net, reshape(A4_10_tmp,128,[])');
    [yy,xx]=find(Results_A4_10_tmp);
%     figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases])
    SOM_A4_10(:,zodzio_nr)=yy;
end
%% Zodyno atnaujinimo sprendimo patikra
% select 4 references
ref_score = zeros(10, zodziu_sk);
ref_max_val = zeros(1,zodziu_sk); ref_max_idx = zeros(1,zodziu_sk);
ref_max_val_pattern = zeros(1,zodziu_sk); ref_max_idx_pattern = zeros(1,zodziu_sk);
ref_max_val_int = zeros(1,zodziu_sk); ref_max_idx_int = zeros(1,zodziu_sk);

% eksperimento tasai
% % load Po_SOM_20140126.mat
SOM_A4_5 = SOM_A4_9;
skaitliukas1=0;
skaitliukas2=0;
for ref_idx = 1:zodziu_sk
    ref_score(1,ref_idx) = nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_5(:,ref_idx))'); 
    ref_score(2,ref_idx) = nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_5(:,ref_idx))'); 
    ref_score(3,ref_idx) = nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_5(:,ref_idx))'); 
    ref_score(4,ref_idx) = nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_5(:,ref_idx))'); 
    ref_score(5,ref_idx) = nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_2(:,ref_idx))'); 
    ref_score(6,ref_idx) = nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_3(:,ref_idx))'); 
    ref_score(7,ref_idx) = nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_4(:,ref_idx))'); 
    ref_score(8,ref_idx) = nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_3(:,ref_idx))'); 
    ref_score(9,ref_idx) = nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_4(:,ref_idx))'); 
    ref_score(10,ref_idx) = nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_4(:,ref_idx))'); 
    [max_val, max_idx] = max(ref_score(:,ref_idx));
    if length(max_idx)>1
        disp('Alarmas!')
    end
    ref_max_val(ref_idx) = max_val(1);
    ref_max_idx(ref_idx) = max_idx(1);

    [max_val_pattern, max_idx_pattern] = max(ref_score(1:4,ref_idx));
    ref_max_val_pattern(ref_idx) = max_val_pattern(1);
    ref_max_idx_pattern(ref_idx) = max_idx_pattern(1);
    
    [max_val_int, max_idx_int] = max(ref_score(5:end,ref_idx));
    ref_max_val_int(ref_idx) = max_val_int(1);
    ref_max_idx_int(ref_idx) = max_idx_int(1);
    if ref_max_val_pattern(ref_idx)<ref_max_val_int(ref_idx)
        disp(['Potencialus zodis atmetimui yra: ', num2str(ref_max_idx_int(ref_idx)+4)])
        skaitliukas1=skaitliukas1+1;
        % tikrinama antra salyga
        
        % potencialaus atmetimui zodzio palyginimas su kitais zodyno
        % zodziais (align two sequences using Needleman-Wunsch algorithm)
        switch ref_max_idx_int(ref_idx)
            case 1
                sum_1 = sum([ref_score(5,ref_idx), ref_score(6,ref_idx), ref_score(7,ref_idx)]);
                sum_2 = sum([ref_score(5,ref_idx), ref_score(8,ref_idx), ref_score(9,ref_idx)]);
                if sum_1 > sum_2
                    voc_score_1 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
                            voc_score_1 = [voc_score_1; nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_1) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 1')
                        A4_1(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                   else 
                    end
                else
                    voc_score_2 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
    %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_2) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 2')
                        A4_2(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
            case 2
                sum_1 = sum([ref_score(5,ref_idx), ref_score(6,ref_idx), ref_score(7,ref_idx)]);
                sum_3 = sum([ref_score(6,ref_idx), ref_score(8,ref_idx), ref_score(10,ref_idx)]);
                if sum_1 > sum_3
                    voc_score_1 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
                            voc_score_1 = [voc_score_1; nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_1) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 1')
                        A4_1(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                else
                    voc_score_3 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_3 = [voc_score_3; nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_3) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 3')
                        A4_1(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
            case 3
                sum_1 = sum([ref_score(5,ref_idx), ref_score(6,ref_idx), ref_score(7,ref_idx)]);
                sum_4 = sum([ref_score(7,ref_idx), ref_score(9,ref_idx), ref_score(10,ref_idx)]);                
                if sum_1 > sum_4
                    voc_score_1 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
                            voc_score_1 = [voc_score_1; nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_1(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_1) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 1')
                        A4_1(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                else
                    voc_score_4 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score4 = [voc_score_4; nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_4) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 4')
                        A4_4(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
            case 4
                sum_2 = sum([ref_score(5,ref_idx), ref_score(8,ref_idx), ref_score(9,ref_idx)]);   
                sum_3 = sum([ref_score(6,ref_idx), ref_score(8,ref_idx), ref_score(10,ref_idx)]);
                if sum_2 > sum_3
                    voc_score_2 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
    %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else 
                        end
                    end
                    if max(voc_score_2) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 2')
                        A4_2(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                else
                    voc_score_3 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_3 = [voc_score_3; nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_3) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 3')
                        A4_3(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
            case 5
                sum_2 = sum([ref_score(5,ref_idx), ref_score(8,ref_idx), ref_score(9,ref_idx)]);
                sum_4 = sum([ref_score(7,ref_idx), ref_score(9,ref_idx), ref_score(10,ref_idx)]);
                if sum_2 > sum_4
                    voc_score_2 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
    %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_2) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 2')
                        A4_2(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                else
                    voc_score_4 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score4 = [voc_score_4; nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_4) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 4')
                        A4_4(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
            case 6 
                sum_3 = sum([ref_score(6,ref_idx), ref_score(8,ref_idx), ref_score(10,ref_idx)]);
                sum_4 = sum([ref_score(7,ref_idx), ref_score(9,ref_idx), ref_score(10,ref_idx)]);
                if sum_3 > sum_4
                    voc_score_3 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score_3 = [voc_score_3; nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_3(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_3) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 3')
                        A4_3(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                else
                    voc_score_4 = [];
                    voc_score_5 = [];
                    for voc_idx = 1:zodziu_sk 
                        if voc_idx ~= ref_idx
        %                     voc_score_2 = [voc_score_2; nwalign(int2aa(SOM_A4_2(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))')]; 
                            voc_score4 = [voc_score_4; nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_4(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')]; 
                            voc_score_5 = [voc_score_5; nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_1(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_2(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_3(:,voc_idx))');...
                                nwalign(int2aa(SOM_A4_5(:,ref_idx))',int2aa(SOM_A4_4(:,voc_idx))')];
                        else
                        end
                    end
                    if max(voc_score_4) > max(voc_score_5)
                        disp('Laikas atnaujinti pakeiciant 4')
                        A4_4(ref_idx,:,:)=A4_5(ref_idx,:,:);
                        skaitliukas2=skaitliukas2+1;
                    else 
                    end
                end
        end
    end
end

    %% DTW atpazinimas
% %     load SOM_10_istarimu_5_naudojami.mat
%     load Ten_closest.mat
% % bendras_int_error = [];
% % for j=1:zodziu_sk
% % %     error = zeros(zodziu_sk,zodziu_sk); word_nr=zeros(1,zodziu_sk); word_nr2=zeros(1,zodziu_sk);
% %     % paieðka kokioje klasëje yra þodis, jei þinome þdþio nr
% % % %     [eil, stulp] = find(Tikslo_mat(:,1)==j);
% % % %     eil;
% % % %     klase = Tikslo_mat(eil,2);
% %      % paieðka kokioje klasëje yra  þodis, jei neþinome þdþio nr
% %      klase = find(sim(net_recognition, SOM_A4(:,j)));
% %     [klases_zodziu_ind, stulp2] = find(Tikslo_mat(:,2)==klase);
% %     klases_zodziai = Tikslo_mat(klases_zodziu_ind,1);
% %     zodziu_sk_klaseje = size(klases_zodziai,1);
% %     for i = 1:zodziu_sk_klaseje% 2:2 % 1:128
% %         serieA = abs(A2(j,1:A2s(j),:)); serieB = abs(A4(klases_zodziai(i),1:A4s(klases_zodziai(i)),:));
% %         lenA = size(serieA);     lenB = size(serieB);
% %         [error(klases_zodziai(i),j),int_error] = dtwV21(serieA,serieB,lenA(2),lenB(2),eile);
% %         bendras_int_error{j,klases_zodziai(i)} = {int_error};
% %     end
% %     error_tmp=error(:,j);
% %     error_tmp(error_tmp==0)=[];
% %     [min_val tmp_value] = min(error_tmp);
% %     word_nr(j)=find(error(:,j)==min_val);
% %     if word_nr(j) == j
% %         cnt_som = cnt_som + 1;
% %         disp([num2str(j),' teisingai'])
% %     else
% %         err(1,k)=j;
% %         err(2,k)=word_nr(j);
% %         k=k+1;
% %         disp([num2str(j),' klaidingai, ',num2str(word_nr(j))])
% %     end
% %     rate_som=100*cnt_som/j;
% % end
% % disp(['Tikslumas ',num2str(rate_som),' %'])
    %% DTW atpazinimas
%     load SOM_10_istarimu_5_naudojami.mat
% % load Po_SOM_20140126.mat
% load pradinis_rezultatas.mat
bendras_rezult = zeros(4,100);
bendras_int_error = [];
bendras_error_visi = [];
matrica = zeros(4,zodziu_sk);
skaitliukas_klaidingai_1 = 0;
skaitliukas_klaidingai_2 = 0;
skaitliukas_klaidingai_3 = 0;
skaitliukas_klaidingai_4 = 0;
for j=1:zodziu_sk
    bendras_error = [];
    for i = 1:zodziu_sk% 2:2 % 1:128
        serieA = abs(A4_10(j,1:A4s(j),:)); 
        serieB_1 = abs(A4_1(i,1:A4s(i),:));
        serieB_2 = abs(A4_2(i,1:A4s(i),:));
        serieB_3 = abs(A4_3(i,1:A4s(i),:));
        serieB_4 = abs(A4_4(i,1:A4s(i),:));
        lenA = size(serieA);     lenB = size(serieB_1);
        [error(i,j),int_error] = dtwV21(serieA,serieB_1,lenA(2),lenB(2),eile);
%         length(int_error)
        bendras_int_error{j,i} = {int_error};
        bendras_error = [bendras_error;error(i,j)];
        [error(i,j),int_error] = dtwV21(serieA,serieB_2,lenA(2),lenB(2),eile);
        bendras_error = [bendras_error;error(i,j)];
        [error(i,j),int_error] = dtwV21(serieA,serieB_3,lenA(2),lenB(2),eile);
        bendras_error = [bendras_error;error(i,j)];
        [error(i,j),int_error] = dtwV21(serieA,serieB_4,lenA(2),lenB(2),eile);
        bendras_error = [bendras_error;error(i,j)];
% %         disp([num2str(j),' ',num2str(i),' ',num2str(rate),'%']);% commented by AS 2013 03 27
    end
%     [min_val word_nr(j)] = min(error(:,j));
%     if word_nr(j) == j
%         cnt = cnt + 1;
%     else
        bendras_error_tmp = bendras_error;
        [min_val word_idx] = min(bendras_error_tmp);
        matrica(1,j) = word_idx;%-4*(j-1);
        bendras_error_tmp(bendras_error_tmp==min(bendras_error_tmp))=1e+09;
        [min_val word_idx] = min(bendras_error_tmp);
        matrica(2,j) = word_idx;%-4*(j-1);
        bendras_error_tmp(bendras_error_tmp==min(bendras_error_tmp))=1e+09;
        [min_val word_idx] = min(bendras_error_tmp);
        matrica(3,j) = word_idx;%-4*(j-1);
        bendras_error_tmp(bendras_error_tmp==min(bendras_error_tmp))=1e+09;
        [min_val word_idx] = min(bendras_error_tmp);
        matrica(4,j) = word_idx;%-4*(j-1);
%         err(1,k)=j;
%         err(2,k)=word_nr(j);
%         k=k+1;
        disp(['Turi buti intervale nuo: ',num2str((j-1)*4+1),' iki ', num2str((j-1)*4+4),' Atpazinti kaip: ', num2str(matrica(1,j)),' ',num2str(matrica(2,j)), ' ',num2str(matrica(3,j)),' ',num2str(matrica(4,j))])
%         disp(['Neteisingai atpaþintas ', num2str(j),' þodis, priskiriant jam ', num2str(word_nr(j)), ' þodá'])
        if (matrica(1,j))>=((j-1)*4+1) && (matrica(1,j))<=((j-1)*4+4) 
        else
            skaitliukas_klaidingai_1=skaitliukas_klaidingai_1+1;
        end
            
        if (matrica(4,j))>=((j-1)*4+1) && (matrica(4,j))<=((j-1)*4+4) 
        else
            skaitliukas_klaidingai_4=skaitliukas_klaidingai_4+1;
        end
        
        if (matrica(2,j))>=((j-1)*4+1) && (matrica(2,j))<=((j-1)*4+4) 
        else
            skaitliukas_klaidingai_2=skaitliukas_klaidingai_2+1;
        end
        if (matrica(3,j))>=((j-1)*4+1) && (matrica(3,j))<=((j-1)*4+4) 
        else
            skaitliukas_klaidingai_3=skaitliukas_klaidingai_3+1;
        end
        bendras_rezult (:,j)=[matrica(1,j);matrica(2,j);matrica(3,j);matrica(4,j)];
        
        %     end
%     [min_val word_nr2(j)] = min(error_norm(:,j));
%     if word_nr2(j) == j
%         cnt2 = cnt2 + 1;
%     else
%         err2(1,k)=j;
%         err2(2,k)=word_nr2(j);
%         k2=k2+1;
%     end
% % %     rate=100*cnt/j;
%     rate2=100*cnt2/j;
    bendras_error_visi = [bendras_error_visi, bendras_error];
end
% rate
disp(['Tikslumas ',num2str(rate),' %'])
