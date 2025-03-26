function [pozymis pozymis0] = LPC_LPCC_v3(signalas, eile)
%
% Isskiriami LPC pozymiai. Naudojamas Levinsono-Durbino algoritmas
%
% Skaiciuojam signalo autokoreliacija
r = xcorr(signalas);
r(1:(find(r==max(r))-1)) = []; % Pasaliname seka iki max
% for i = 1:13
%     r(i) = r(i)/r(1);
% end
% Skaiciuojam LPC koeficientus
temp = 0;
ki = zeros(1, eile); % Dalines koreliacijos koeficientai
lpcKoef = zeros(1, eile); % LPC koeficientai
a = lpcKoef; % Laikinas
lpcKoefTemp2 = lpcKoef; % Laikinas

E = r(1);
ki(1) = r(1+1) / E; % I-asis dalines koreliacijos koef.
a(1) = ki(1); % I-asis LPC koef.
E = (1 - ki(1)^2) * E; % Energijos reiksme
for i = 2:eile 
    temp = 0;
    for j = 1:(i-1) % Skaiciuojam dalines koreliacijos koef.
        temp = temp + a(j)*r(i-j+1);
    end;
    ki(i) = (r(i+1) - temp) / E; % dalines kor. koef.
%     if(abs(ki(i))>=1)
%         disp('asdfasdfasdf'); [y] = wavread(pavad); 
%     end
    lpcKoefTemp2(i) = ki(i); % i-asis LPC koeficientas
    for j = 1:(i-1)
        lpcKoefTemp2(j) = a(j) - ki(i) * a(i-j);
    end;
    a = lpcKoefTemp2;
    E = (1 - ki(i)^2) * E; % Nauja energijos reiksme
end;
% pozymis = (lpcKoefTemp2);
a = lpcKoefTemp2;
% pozymis = fix(lpcKoefTemp2*1000);
% figure(4); bar(pozymis); grid on; axis([0 eile+1 -3 3]);
c = zeros(1,12);
c = a;
for i=2:12
    for j=1:i-1
        c(i) = c(i) + (j)*c(j)*a(i-j)/i;
    end
end
pozymis = (c); pozymis0 = (a);
figure(2); bar(pozymis); grid on; axis([0 eile+1 -0.5 3]);