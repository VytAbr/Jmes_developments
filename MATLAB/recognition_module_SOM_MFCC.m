clc; 
clear all; 
close all;
zodziu_sk = 100; s=zeros(128,256); eile = 12;
A2 = zeros(zodziu_sk,128,eile);  C2 = zeros(zodziu_sk,128,eile); 
A2s = zeros(1,zodziu_sk);
A4 = zeros(zodziu_sk,128,eile);  C4 = zeros(zodziu_sk,128,eile); 
A4s = zeros(1,zodziu_sk);
for j=1:zodziu_sk % 1:100
%     [data_y4 length_y] = pavad_normalization_0215(j,'/M030',3);
    [data_y4 length_y] = pavad_normalization_0405(j,'/V030',5);
    [data_y2 length_y] = pavad_normalization_0405(j,'/V030',7);
    len_y = floor(length_y);
    for i = 1:128 % 2:2 % 1:128
% % % % % % % % -----------------------------------------
%         s=data_y4((1+(i-1)*256):(i*256));
        s = fix((data_y4((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
%         s_tr = fft(data_y4(1:256));
        s_tr = 0;

%         figure(2);plot(s); grid on; axis([0 256 -128 128]);
% %         [x A4(j,i,:) ] = LPC_LPCC_v3(s,eile); a4(:,:)=A4(j,i,:);
%           A4(j,i,:) = cepstrum_20130405(s,eile,s_tr); a4(:,:)=A4(j,i,:);
% % %         figure(4); 
% % %         bar(a4,'r'); grid on; axis([0 eile+1 -4 4]); title('LPC');
%         figure(6); 
%         bar(x,'r'); grid on; axis([0 eile+1 -4 4]); title('x');
%         figure(5); 
%         bar(c4,'r'); grid on; axis([0 eile+1 -3 1000]); title('Cepstrum');
        A4(j,i,:) = mfcc(s,128); a4(:,:)=A4(j,i,:);
        A4s(j) = len_y;
% % % % % % % % -----------------------------------------
%         t=data_y6((1+(i-1)*256):(i*256));
        t = fix((data_y2((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
        t_tr = 0;

%         figure(2);plot(t); grid on; axis([0 256 -128 128]);
% %         [y A2(j,i,:) ] = LPC_LPCC_v3(t,eile);   a2(:,:)=A2(j,i,:);
%         A2(j,i,:) = cepstrum_20130405(t,eile,t_tr);       a2(:,:)=A2(j,i,:);
% % %         figure(4); 
% % %         bar(a2,'g'); grid on; axis([0 eile+1 -4 4]); title('LPC');
%         figure(6); 
%         bar(y,'g'); grid on; axis([0 eile+1 -4 4]); title('y');
%         figure(5); 
%         bar(c2,'g'); grid on; axis([0 eile+1 -3 1000]); title('Cepstrum');
        A2(j,i,:) = mfcc(t,128);       a2(:,:)=A2(j,i,:);
        A2s(j) = len_y;
%         figure(6);
%         bar(abs(c4-c2),'g'); grid on; axis([0 eile+1 0 200]); title('Skirtumas Cepstrum');
%         figure(7);
%         bar(abs(a4-a2),'g'); grid on; axis([0 eile+1 0 2]); title('Skirtumas LPC');
% % % % % % % % -----------------------------------------
% %         C1(j,i,:) = cepstrum(s(i,:));
%         disp(i);
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
net = train(net,reshape(A2,12800,[])');
% % Results_A2 = sim(net, reshape(A2,12800,[])');
% % Results_A4 = sim(net, reshape(A4,12800,[])');
% % [yy,xx]=find(Results_A2);
close all
SOM_A2=zeros(128,100);
SOM_A4=zeros(128,100);
for zodzio_nr = 1%:zodziu_sk
    A2_tmp = A2(zodzio_nr,:,:);
    Results_A2_tmp = sim(net, reshape(A2_tmp,128,[])');
    [yy,xx]=find(Results_A2_tmp);
    figure(zodzio_nr),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases(1)]), grid on
    SOM_A2(:,zodzio_nr)=yy;
end
for zodzio_nr = 1%:zodziu_sk
    A4_tmp = A4(zodzio_nr,:,:);
    Results_A4_tmp = sim(net, reshape(A4_tmp,128,[])');
    [yy,xx]=find(Results_A4_tmp);
    figure(zodzio_nr+6),plot(xx(1:128),yy(1:128),'k*-'),axis([0 128 0 klases(1)]), grid on
    SOM_A4(:,zodzio_nr)=yy;
end
%% Klasifikavimas pagal SOM isejima
% 5 triukasma panasios pradzios naikinimas - 5 klasei priskirtu pozymiu
% eliminavimas
pirma_kl_A2 = [];
antra_kl_A2 = [];
trecia_kl_A2 = [];
ketvirta_kl_A2 = [];
sesta_kl_A2 = [];

pirma_kl_A4 = [];
antra_kl_A4 = [];
trecia_kl_A4 = [];
ketvirta_kl_A4 = [];
sesta_kl_A4 = [];

zodziai_A2=[];
zodziai_A4=[];
sekos_A2=[];
sekos_A4=[];
for indeksas=1:zodziu_sk
    zodis=[];
    r=SOM_A2(:,indeksas)';
    ii=1;
    while ii<(numel(r)-2)
        if r(ii)==r(ii+1) && r(ii)==r(ii+2) && r(ii)==r(ii+3)
            zodis=[zodis, r(ii), r(ii+1), r(ii+2), r(ii+3)];
            ii=ii+4;
        else
            ii=ii+1;
        end
    end
    zodziai_A2{indeksas}=zodis;

    
    %% Skirstymas klasemis
    % reikia ieskoti bent keturiu einanciu vienas po kito tos pacios klases
    % elementu
    seka_A2=[];
    while isempty(zodis)==0
        switch zodis(1)
            case 1
                seka_A2=[seka_A2,1];
                while isempty(zodis)==0 && zodis(1)==1
                    zodis(1)=[];
                end
%                 zodis(1)=[];
            case 2
                seka_A2=[seka_A2,2];
                while isempty(zodis)==0 && zodis(1)==2
                    zodis(1)=[];
                end
%                 zodis(1)=[];
            case 3
                seka_A2=[seka_A2,3];
                while isempty(zodis)==0 && zodis(1)==3
                    zodis(1)=[];
                end
%                 zodis(1)=[];
            case 4
                seka_A2=[seka_A2,4];
                while isempty(zodis)==0 && zodis(1)==4
                    zodis(1)=[];
                end
%                 zodis(1)=[];
            case 5
                seka_A2=[seka_A2,5];
                
                while isempty(zodis)==0 && zodis(1)==5
                    zodis(1)=[];
                end
%                 zodis(1)=[];
            case 6
                seka_A2=[seka_A2,6];
                while isempty(zodis)==0 && zodis(1)==6
                    zodis(1)=[];
                end
%                 zodis(1)=[];
        end
    end
    sekos_A2{indeksas}=seka_A2;
% %     while r(end)==5
% %         r(end)=[];
% %     end    
% % %     if zodis(1)==1 || zodis(2)==1
% % %         pirma_kl_A2 = [pirma_kl_A2; indeksas];
% % %     else if zodis(1)==2 || zodis(2)==2
% % %             antra_kl_A2 = [antra_kl_A2; indeksas];
% % %         else if zodis(1)==3 || zodis(2)==3
% % %                 trecia_kl_A2 = [trecia_kl_A2; indeksas];
% % %             else if zodis(1)==4 || zodis(2)==4
% % %                     ketvirta_kl_A2 = [ketvirta_kl_A2; indeksas];
% % %                 else
% % %                     sesta_kl_A2 = [sesta_kl_A2; indeksas];
% % %                 end
% % %             end
% % %         end
% % %     end
    
    zodis=[];
    r=SOM_A4(:,indeksas)';
    ii=1;
    while ii<(numel(r)-2)
        if r(ii)==r(ii+1) && r(ii)==r(ii+2) && r(ii)==r(ii+3)
            zodis=[zodis, r(ii), r(ii+1), r(ii+2), r(ii+3)];
            ii=ii+4;
        else
            ii=ii+1;
        end
    end
    zodziai_A4{indeksas}=zodis;

seka_A4=[];
    while isempty(zodis)==0
        switch zodis(1)
            case 1
                seka_A4=[seka_A4,1];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==1
                    zodis(1)=[];
                end
            case 2
                seka_A4=[seka_A4,2];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==2
                    zodis(1)=[];
                end
            case 3
                seka_A4=[seka_A4,3];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==3
                    zodis(1)=[];
                end
            case 4
                seka_A4=[seka_A4,4];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==4
                    zodis(1)=[];
                end
            case 5
                seka_A4=[seka_A4,5];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==5 
                    zodis(1)=[];
                end
            case 6
                seka_A4=[seka_A4,6];
%                 zodis(1)=[];
                while isempty(zodis)==0 && zodis(1)==6
                    zodis(1)=[];
                end
        end
    end
sekos_A4{indeksas}=seka_A4;
% % %     while zodis(1)==5
% % %         zodis(1)=[];
% % %     end
% %     while r(end)==5
% %         r(end)=[];
% %     end    
% % %     if zodis(1)==1 || zodis(2)==1
% % %         pirma_kl_A4 = [pirma_kl_A4; indeksas];
% % %     else if zodis(1)==2 || zodis(2)==2
% % %             antra_kl_A4 = [antra_kl_A4; indeksas];
% % %         else if zodis(1)==3 || zodis(2)==3
% % %                 trecia_kl_A4 = [trecia_kl_A4; indeksas];
% % %             else if zodis(1)==4 || zodis(2)==4
% % %                     ketvirta_kl_A4 = [ketvirta_kl_A4; indeksas];
% % %                 else
% % %                     sesta_kl_A4 = [sesta_kl_A4; indeksas];
% % %                 end
% % %             end
% % %         end
% % %     end
    
end
%% Globally align two sequences using Needleman-Wunsch algorithm
% % % A4_sequences = [];
% % % % score = [];
% % % for indeksas_align = 1:zodziu_sk
% % %     score = [];
% % %     r=SOM_A2(:,indeksas_align)';
% % %     A2_test=int2aa(r);% tikslumas 68%
% % % %     A2_test=int2aa(zodziai_A2{indeksas_align});
% % %     for indeksas_word = 1:zodziu_sk
% % %         A4_sequences{indeksas_word} = int2aa(SOM_A4(:,indeksas_word)');
% % % %         A4_sequences{indeksas_word} = int2aa(zodziai_A4{indeksas_word});
% % %         score = [score;nwalign(A4_sequences{indeksas_word},A2_test)];
% % %     end
% % %     [score_a,score_b]=max(score);
% % %     disp(score_b)
% % % end

A4_sequences = [];
% score = [];
cnt3=0;
cnt4=0;
for indeksas_align = 1:zodziu_sk
    score = [];
    Score = [];
    A2_test=int2aa(sekos_A2{indeksas_align});% tikslumas 62 %
% %     r=SOM_A2(:,indeksas_align)';
% %     A2_test=int2aa(r);% tikslumas 68%
%     A2_test=int2aa(zodziai_A2{indeksas_align});
    for indeksas_word = 1:zodziu_sk
        A4_sequences{indeksas_word} = int2aa(sekos_A4{indeksas_word});
% %         A4_sequences{indeksas_word} = int2aa(SOM_A4(:,indeksas_word)');
%         A4_sequences{indeksas_word} = int2aa(zodziai_A4{indeksas_word});
        score = [score;nwalign(A4_sequences{indeksas_word},A2_test,'GapOpen',5)];
%         score = [score;swalign(A4_sequences{indeksas_word},A2_test,'GapOpen',5)];
%         score{indeksas_word} = localalign(A4_sequences{indeksas_word},A2_test,'GapOpen',4);
%         Score=[Score,score{1,indeksas_word}.Score];
    end
    [score_a,score_b]=max(score);
%     [score_a,score_b]=max(Score);
    [find_a,find_b]=find(score>7);
    if score_b==indeksas_align
        cnt3=cnt3+1;
    else
    end
    [find_aa,find_bb]=find(score>9);
    disp(numel(find_aa))
    if find(find_aa==indeksas_align)
        cnt4=cnt4+1;
    else
    end
%     disp(score_b)
end
disp(['Tikslumas ', num2str(cnt3),' %'])
disp(['Tikslumas ', num2str(cnt4),' %'])
%% Atvaizdavimas
% % % for indeksas_plot=1:8:zodziu_sk
% % %     figure,
% % %     subplot(441), plot(1:numel(zodziai_A2{indeksas_plot}), zodziai_A2{indeksas_plot}),axis([0 110 0 9])
% % %     subplot(445), plot(1:numel(zodziai_A4{indeksas_plot}), zodziai_A4{indeksas_plot}),axis([0 110 0 9])
% % %     subplot(442), plot(1:numel(zodziai_A2{indeksas_plot+1}), zodziai_A2{indeksas_plot+1}),axis([0 110 0 9])
% % %     subplot(446), plot(1:numel(zodziai_A4{indeksas_plot+1}), zodziai_A4{indeksas_plot+1}),axis([0 110 0 9])
% % %     subplot(443), plot(1:numel(zodziai_A2{indeksas_plot+2}), zodziai_A2{indeksas_plot+2}),axis([0 110 0 9])
% % %     subplot(447), plot(1:numel(zodziai_A4{indeksas_plot+2}), zodziai_A4{indeksas_plot+2}),axis([0 110 0 9])
% % %     subplot(444), plot(1:numel(zodziai_A2{indeksas_plot+3}), zodziai_A2{indeksas_plot+3}),axis([0 110 0 9])
% % %     subplot(448), plot(1:numel(zodziai_A4{indeksas_plot+3}), zodziai_A4{indeksas_plot+3}),axis([0 110 0 9])
% % %     subplot(449), plot(1:numel(zodziai_A2{indeksas_plot+4}), zodziai_A2{indeksas_plot+4}),axis([0 110 0 9])
% % %     subplot(4,4,13), plot(1:numel(zodziai_A4{indeksas_plot+4}), zodziai_A4{indeksas_plot+4}),axis([0 110 0 9])
% % %     subplot(4,4,10), plot(1:numel(zodziai_A2{indeksas_plot+5}), zodziai_A2{indeksas_plot+5}),axis([0 110 0 9])
% % %     subplot(4,4,14), plot(1:numel(zodziai_A4{indeksas_plot+5}), zodziai_A4{indeksas_plot+5}),axis([0 110 0 9])
% % %     subplot(4,4,11), plot(1:numel(zodziai_A2{indeksas_plot+6}), zodziai_A2{indeksas_plot+6}),axis([0 110 0 9])
% % %     subplot(4,4,15), plot(1:numel(zodziai_A4{indeksas_plot+6}), zodziai_A4{indeksas_plot+6}),axis([0 110 0 9])
% % %     subplot(4,4,12), plot(1:numel(zodziai_A2{indeksas_plot+7}), zodziai_A2{indeksas_plot+7}),axis([0 110 0 9])
% % %     subplot(4,4,16), plot(1:numel(zodziai_A4{indeksas_plot+7}), zodziai_A4{indeksas_plot+7}),axis([0 110 0 9])
% % % end
% % for indeksas=1:zodziu_sk
% %     r=SOM_A4(:,indeksas)';
% %     while r(1)==5
% %         r(1)=[];
% %     end
% %     while r(end)==5
% %         r(end)=[];
% %     end
% %     if r(1)==1 || r(2)==1
% %         if find(pirma_kl==indeksas);
% %             disp([num2str(indeksas), ' suskirstytas teisingai']);
% %         else
% %             disp('Klaida')
% %         end
% %     else if r(1)==2 || r(2)==2
% % %             antra_kl = [antra_kl; indeksas];
% %             if find(antra_kl==indeksas);
% %                 disp([num2str(indeksas), ' suskirstytas teisingai']);
% %             else
% %                 disp('Klaida')
% %             end
% %         else if r(1)==3 || r(2)==3
% % %                 trecia_kl = [trecia_kl; indeksas];
% %                 if find(trecia_kl==indeksas);
% %                     disp([num2str(indeksas), ' suskirstytas teisingai']);
% %                 else
% %                     disp('Klaida')
% %                 end
% %             else if r(1)==4 || r(2)==4
% % %                     ketvirta_kl = [ketvirta_kl; indeksas];
% %                     if find(ketvirta_kl==indeksas);
% %                         disp([num2str(indeksas), ' suskirstytas teisingai']);
% %                     else
% %                         disp('Klaida')
% %                     end
% %                 else
% % %                     sesta_kl = [sesta_kl; indeksas];
% %                     if find(sesta_kl==indeksas);
% %                         disp([num2str(indeksas), ' suskirstytas teisingai']);
% %                     else
% %                         disp('Klaida')
% %                     end
% %                 end
% %             end
% %         end
% %     end
% %     
% % end

%% SOM atpaþinimas
% % % net_recognition = selforgmap([10 1]);
% % % net_recognition = train(net_recognition,SOM_A2);
% % % % Results_A2_rec = sim(net_recognition, SOM_A2);
% % % % [yy,xx]=find(Results_A2_rec);
% % % Results_A4_rec = sim(net_recognition, SOM_A4);
% % % [yy,xx]=find(Results_A4_rec);
% % % Tikslo_mat = sortrows([xx,yy],2);
%% DTW atpaþinimas ið SOM
% % % Dist=zeros(zodziu_sk,1);
% % % skaitiklis=0;
% % % skait=0;
% % % Nr=zeros(zodziu_sk,1);
% % % Ten_closest = zeros(zodziu_sk,10);
% % % for indeksas=1:zodziu_sk
% % % %     r=SOM_A2(:,indeksas)';
% % %     r=zodziai_A2{indeksas};
% % % %     while r(1)==5
% % % %         r(1)=[];
% % % %     end
% % % %     while r(end)==5
% % % %         r(end)=[];
% % % %     end
% % %     for zodyno_indeksas=1:zodziu_sk
% % % %         t=SOM_A4(:,zodyno_indeksas)';
% % %         t=zodziai_A4{zodyno_indeksas};
% % % %         while t(1)==5
% % % %             t(1)=[];
% % % %         end
% % % %         while t(end)==5
% % % %             t(end)=[];
% % % %         end
% % % %         [Dist(zodyno_indeksas),D,k,w]=dtw_TF(medfilt1(t,3),medfilt1(r,3));%filtravimas
% % %         [Dist(zodyno_indeksas),D,k,w]=dtw_TF(t,r);
% % %     end
% % %     [b,a]=min(Dist);
% % %     if a==indeksas
% % %         disp('Correct')
% % %         Dist2=Dist;
% % %         Dist2(a)=max(Dist);
% % %         [b2,a2]=min(Dist2);
% % %         Dist3=Dist2;
% % %         Dist3(a2)=max(Dist2);
% % %         [b3,a3]=min(Dist3);
% % %         Dist4=Dist3;
% % %         Dist4(a3)=max(Dist3);
% % %         [b4,a4]=min(Dist4);
% % %         Dist5=Dist4;
% % %         Dist5(a4)=max(Dist4);
% % %         [b5,a5]=min(Dist5);
% % %         Dist6=Dist5;
% % %         Dist6(a5)=max(Dist5);
% % %         [b6,a6]=min(Dist6);
% % %         
% % %         Dist7=Dist6;
% % %         Dist7(a6)=max(Dist6);
% % %         [b7,a7]=min(Dist7);
% % %         
% % %         Dist8=Dist7;
% % %         Dist8(a7)=max(Dist7);
% % %         [b8,a8]=min(Dist8);
% % % 
% % %         Dist9=Dist8;
% % %         Dist9(a8)=max(Dist8);
% % %         [b9,a9]=min(Dist9);
% % % 
% % %         Dist10=Dist9;
% % %         Dist10(a9)=max(Dist9);
% % %         [b10,a10]=min(Dist10);
% % %         Ten_closest(indeksas,:) = [a,a2,a3,a4,a5,a6,a7,a8,a9,a10];
% % %     else
% % %         Dist2=Dist;
% % %         Dist2(a)=max(Dist);
% % %         [b2,a2]=min(Dist2);
% % %         Dist3=Dist2;
% % %         Dist3(a2)=max(Dist2);
% % %         [b3,a3]=min(Dist3);
% % %         Dist4=Dist3;
% % %         Dist4(a3)=max(Dist3);
% % %         [b4,a4]=min(Dist4);
% % %         Dist5=Dist4;
% % %         Dist5(a4)=max(Dist4);
% % %         [b5,a5]=min(Dist5);
% % %         Dist6=Dist5;
% % %         Dist6(a5)=max(Dist5);
% % %         [b6,a6]=min(Dist6);
% % %         
% % %         Dist7=Dist6;
% % %         Dist7(a6)=max(Dist6);
% % %         [b7,a7]=min(Dist7);
% % %         
% % %         Dist8=Dist7;
% % %         Dist8(a7)=max(Dist7);
% % %         [b8,a8]=min(Dist8);
% % % 
% % %         Dist9=Dist8;
% % %         Dist9(a8)=max(Dist8);
% % %         [b9,a9]=min(Dist9);
% % % 
% % %         Dist10=Dist9;
% % %         Dist10(a9)=max(Dist9);
% % %         [b10,a10]=min(Dist10);
% % %         Ten_closest(indeksas,:) = [a,a2,a3,a4,a5,a6,a7,a8,a9,a10];
% % %         disp(['error at ', num2str(indeksas), ' found as ', num2str(a),...
% % %             ' Other: ', num2str(a2),' ', num2str(a3), ' ',num2str(a4),...
% % %             ' ',num2str(a5), ' ',num2str(a6), ' ',num2str(a7), ' ',num2str(a8), ' ',num2str(a9), ' ',num2str(a10)])
% % %         if a2==indeksas || a3==indeksas || a4==indeksas || a5==indeksas ||...
% % %                 a6==indeksas || a7==indeksas || a8==indeksas || a9==indeksas || a10==indeksas
% % %             skait=skait+1;
% % %             [reiksme,Nr(indeksas)]=find([a,a2,a3,a4,a5,a6,a7,a8,a9,a10]==indeksas);
% % %         else
% % %             
% % %         end
% % %         skaitiklis=skaitiklis+1;
% % %     end
% % % end
% % % tikslumas=100*(zodziu_sk-skaitiklis)/zodziu_sk;
% % % tikslumas2=100*(zodziu_sk-skaitiklis+skait)/zodziu_sk;
% % % disp(['Tikslumas ',num2str(tikslumas),' %'])
% % % disp(['Tikslumas praplëtus ribas ',num2str(tikslumas2),' %'])
% % % max(Nr)
    %% DTW atpazinimas
%     load Ten_closest.mat
% % bendras_int_error = [];
% % for j=1:zodziu_sk
% % %     error = zeros(zodziu_sk,zodziu_sk); word_nr=zeros(1,zodziu_sk); word_nr2=zeros(1,zodziu_sk);
% %     % paieðka kokioje klasëje yra þodis, jei þinome þdþio nr
% % % %     [eil, stulp] = find(Tikslo_mat(:,1)==j);
% % % %     eil;
% % % %     klase = Tikslo_mat(eil,2);
% %      % paieðka kokioje klasëje yra þodis, jei neþinome þdþio nr
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
% % % bendras_int_error = [];
% % % for j=1:zodziu_sk
% % %     for i = 1:zodziu_sk% 2:2 % 1:128
% % %         serieA = abs(A2(j,1:A2s(j),:)); serieB = abs(A4(i,1:A4s(i),:));
% % %         lenA = size(serieA);     lenB = size(serieB);
% % %         [error(i,j),int_error] = dtwV21(serieA,serieB,lenA(2),lenB(2),eile);
% % % %         length(int_error)
% % %         bendras_int_error{j,i} = {int_error};
% % % % %         disp([num2str(j),' ',num2str(i),' ',num2str(rate),'%']);% commented by AS 2013 03 27
% % %     end
% % %     [min_val word_nr(j)] = min(error(:,j));
% % %     if word_nr(j) == j
% % %         cnt = cnt + 1;
% % %     else
% % %         err(1,k)=j;
% % %         err(2,k)=word_nr(j);
% % %         k=k+1;
% % %         disp(['Neteisingai atpaþintas ', num2str(j),' þodis, priskiriant jam ', num2str(word_nr(j)), ' þodá'])
% % %     end
% % % %     [min_val word_nr2(j)] = min(error_norm(:,j));
% % % %     if word_nr2(j) == j
% % % %         cnt2 = cnt2 + 1;
% % % %     else
% % % %         err2(1,k)=j;
% % % %         err2(2,k)=word_nr2(j);
% % % %         k2=k2+1;
% % % %     end
% % %     rate=100*cnt/j;
% % % %     rate2=100*cnt2/j;
% % % end
% % % % rate
% % % disp(['Tikslumas ',num2str(rate),' %'])
%% Rezultato atvaizdavimas
% % zodziu_sk = 15
% for j=1:zodziu_sk
%     for i=1:zodziu_sk
%         if i==word_nr(j)
%             figure(j)
%             hold on
%             plot(cell2mat(bendras_int_error{j,i}),'green')
%             hold off
%         else
%             figure(j)
%             hold on
%             plot(cell2mat(bendras_int_error{j,i}))
%             hold off
%         end
%     end
% end
% % for j=1:zodziu_sk
% %     for i=1:zodziu_sk
% %         if i==word_nr2(j)            
% %             figure(j)
% %             hold on
% %             plot(cell2mat(bendras_int_error{j,i}),'black')
% %             hold off
% %         else
% %             figure(j)
% %             hold on
% %             plot(cell2mat(bendras_int_error{j,i}))
% %             hold off
% %         end
% %         if i==word_nr2(j) & i==j
% %             figure(j)
% %             hold on
% %             plot(cell2mat(bendras_int_error{j,i}),'red')
% %             hold off
% %         else 
% %             if i==j
% %             figure(j)
% %             hold on
% %             plot(cell2mat(bendras_int_error{j,i}),'green')
% %             hold off
% %             else
% %             end
% %         end
% %     end
% % end
