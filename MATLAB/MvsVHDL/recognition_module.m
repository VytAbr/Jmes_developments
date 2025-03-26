clc; 
clear all; 
zodziu_sk = 100; s=zeros(128,256); eile = 12;
A2 = zeros(zodziu_sk,128,eile);  C2 = zeros(zodziu_sk,128,eile); 
A2s = zeros(1,zodziu_sk);
A4 = zeros(zodziu_sk,128,eile);  C4 = zeros(zodziu_sk,128,eile); 
A4s = zeros(1,zodziu_sk);
for j=6:6
%     [data_y4 length_y] = pavad_normalization(j,3);
    [data_y2 length_y] = pavad_normalization(j,1);
    len_y = floor(length_y);
    for i = 1:100 % 2:2 % 1:128
% % % % % % % % -----------------------------------------
% %         s=data_y4((1+(i-1)*256):(i*256));
%         s = ((data_y4((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
% %         figure(2);plot(s); grid on; axis([0 256 -128 128]);
%         [A4(j,i,:) x ] = LPC_LPCC_v3(s,eile); a4(:,:)=A4(j,i,:);
% %           A4(j,i,:) = cepstrum(s,eile); a4(:,:)=A4(j,i,:);
% %         figure(4); 
% %         bar(a4,'r'); grid on; axis([0 eile+1 -4 4]); title('LPC');
% %         figure(6); 
% %         bar(x,'r'); grid on; axis([0 eile+1 -4 4]); title('x');
% %         figure(5); 
% %         bar(c4,'r'); grid on; axis([0 eile+1 -3 1000]); title('Cepstrum');
%         A4s(j) = len_y;
% % % % % % % % -----------------------------------------
%         t=data_y6((1+(i-1)*256):(i*256));
        t = ((data_y2((1+(i-1)*128):(i*128+128))).*hamming(256)'); % .*hann(256)'
%         figure(2);plot(t); grid on; axis([0 256 -128 128]);
        [A2(j,i,:) y ] = LPC_LPCC_v3(t,eile);   a2(:,:)=A2(j,i,:);
%         A2(j,i,:) = cepstrum(t,eile);       a2(:,:)=A2(j,i,:);
%         figure(4); 
%         bar(a2,'g'); grid on; axis([0 eile+1 -4 4]); title('LPC');
%         figure(6); 
%         bar(y,'g'); grid on; axis([0 eile+1 -4 4]); title('y');
%         figure(5); 
%         bar(c2,'g'); grid on; axis([0 eile+1 -3 1000]); title('Cepstrum');
        A2s(j) = len_y;
%         figure(6);
%         bar(abs(c4-c2),'g'); grid on; axis([0 eile+1 0 200]); title('Skirtumas Cepstrum');
%         figure(7);
%         bar(abs(a4-a2),'g'); grid on; axis([0 eile+1 0 2]); title('Skirtumas LPC');
% % % % % % % % -----------------------------------------
% %         C1(j,i,:) = cepstrum(s(i,:));
%         disp(i);
    end
    disp(j);
end
% figure(7);
% plot(t); grid on;

% error = zeros(100,100); word_nr=zeros(1,100); cnt = 0; rate = 0; k = 1;
% % serieA = zeros(128,12); serieB = zeros(128,12);
% for j=1:100
%     for i = 1:100 % 2:2 % 1:128
%         serieA = A2(j,1:A2s(j),:); serieB = A4(i,1:A4s(i),:);
%         lenA = size(serieA);     lenB = size(serieB);
%         error(i,j) = dtw(serieA,serieB,lenA(2),lenB(2),eile);
%         disp([num2str(j),' ',num2str(i),' ',num2str(rate),'%']);
%     end
%     [min_val word_nr(j)] = min(error(:,j));
%     if word_nr(j) == j
%         cnt = cnt + 1;
%     else
%         err(1,k)=j;
%         err(2,k)=word_nr(j);
%         k=k+1;
%     end
%     rate=100*cnt/j;
% end
% rate