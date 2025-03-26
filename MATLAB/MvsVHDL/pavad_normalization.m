function [data_y length_y] = pavad_normalization(i,bank)
% % --------------------------------
pavad3 = '.wav';
pavad1 = cat(2,'V',num2str(bank));
pavad1 = cat(2,pavad1,'_apmokymui/V010');
pavad1 = cat(2,pavad1,num2str(bank));
if i<10
    pavad1 = cat(2,pavad1,'00'); % 'V3_apmokymui/V010300'; 
elseif i<100
    pavad1 = cat(2,pavad1,'0'); % 'V3_apmokymui/V01030';
end
pavad = cat(2,pavad1,num2str(i)); 
pavad = cat(2,pavad,pavad3);
% % --------------------------------
% pavad3 = '.wav';
% pavad1 = cat(2,'15dB/M030',num2str(bank));
% pavad1 = cat(2,pavad1,'/!M030');
% pavad1 = cat(2,pavad1,num2str(bank));
% if i<10
%     pavad1 = cat(2,pavad1,'00'); % 'V3_apmokymui/V010300'; 
% elseif i<100
%     pavad1 = cat(2,pavad1,'0'); % 'V3_apmokymui/V01030';
% end
% pavad = cat(2,pavad1,num2str(i)); 
% pavad = cat(2,pavad,pavad3);
% % ---------------------------------
[y] = wavread(pavad); 
wavplay(y*5,11025); % ------------------- *.wav play
y_max = max(y); y_min = min(y);
d_len = 128*128; % d_len = 16384; 
if abs(y_max) > abs(y_min)
    coef = 128/abs(y_max); % upper_lim = 0.08; 0.7; 
else
    coef = 128/abs(y_min); % upper_lim = 0.08; 0.7; 
end
% coef = 500;% % % % disp(coef);
y_len = length(y);
% -------------------------------------------- tylos nesalinimas
if mod(y_len,2)==1
    o=2; 
else
    o=1; 
end; 
y_centre = round(y_len/2);
y2 = zeros(1,50000)+0.002; y2(25000-y_centre+o:25000+y_centre) = y;
data_y = (coef*y2((25000-d_len/2+1-64):(25000+d_len/2+64))+0);
length_y = 128; data_y_len = length(data_y);
% -------------------------------------------- tylos pasalinimas start
% data_y_tmp = (coef*y-3)';
% I = floor(y_len/128);
% for i = 1:I
%     s=data_y_tmp((1+(i-1)*128):(i*128+128));
%     if sum(abs(s)) > 800
%     	lim_left = i;
%         break;
%     end
% end
% for j = 1:I
%     s=data_y_tmp((1+(I-j)*128):((I-j)*128+128));
%     if sum(abs(s)) > 800
%     	lim_right = j;
%         break;
%     end
% end
% length_y = I - lim_left - lim_right;
% data_y = data_y_tmp((lim_left)*128:(I-lim_right+1)*128-1); data_y_len = length(data_y);
% -------------------------------------------- tylos pasalinimas end
figure(1);plot(data_y); grid on; axis([0 data_y_len -128 128]);