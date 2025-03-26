% clear all; clc;
function error = dtw(serieA,serieB,lenA,lenB,eile)
% lenA = 128;
% lenB = 128;
% serieA = ones(lenA,12)*1; % signalas A
% serieB = ones(lenB,12)*0; % *546; % *273; % signalas B
% % % dtw skaiciuojamas tarp serieA ir serieB signalu
d = zeros(lenB, lenA); g = zeros(lenB, lenA); diff = zeros(1,12);
for b=1:lenB
    for a=1:lenA 
        for c=1:eile
            diff(1,c) = (serieA(1,a,c) - serieB(1,b,c)).^2;
        end
%       d(b,a) = floor(sqrt(sum(diff)));
        d(b,a) = (sqrt(sum(diff)));
    end
end
% d - tai klaidos pavirsius
% g - tai perskaiciuotas klaidos pavirsius, kai sekantis
% g(i,j) elementas sudarytas is d + min(3 kaimyniniai elementai)
g(1,1) = d(1,1);
for i = 2:lenA % uzpildoma eilute 1
    g(1,i) = g(1,i-1)+d(1,i);
end
for j = 2:lenB % uzpildoma stulpelis 1
    g(j,1) = g(j-1,1)+d(j,1);
end
for j = 2:lenB  % uzpildomas likes plotas
    for i = 2:lenA
        g(j,i) = min([g(j-1,i),g(j-1,i-1),g(j,i-1)])+d(j,i);
    end
end
% DTW (kelio, su maziausia kaupiama klaida nuo tasko(128,128)
% iki pradzios tasko(0,0) radimas)
pathArea = zeros(lenA, lenB); 
pathArea(lenA,lenB) = 2;
ax = lenA; bx = lenB;
errorSum = g(bx,ax); 
gmax = max(max(g))*0.75; i = 1;
while (ax > 1 || bx > 1)
        if     ax == 1
            bx = bx - 1; add_error = g(bx,ax);
        elseif bx == 1
            ax = ax - 1; add_error = g(bx,ax);
        elseif (g(bx-1,ax-1) <= g(bx-1,ax) && g(bx-1,ax-1) <= g(bx,ax-1))
            ax = ax - 1; bx = bx - 1; add_error = g(bx,ax);
        elseif (g(bx,ax-1) < g(bx-1,ax-1) && g(bx,ax-1) <= g(bx-1,ax))
            ax = ax - 1; add_error = g(bx,ax);
        else
            bx = bx - 1; add_error = g(bx,ax);
        end
        Gpath(i)=errorSum + g(bx,ax); i=i+1;
        pathArea(bx,ax) = gmax; errorSum = errorSum + add_error; % klaidos kaupimas
end
% disp(['Error sum dec: ', num2str(errorSum)]);
% disp(['Error sum hex: ', num2str(dec2hex(errorSum))]); 
error = errorSum; % ---- atsakymas
% figure(1); subplot(2,2,1); plot(serieA); title('serieA'); grid on;
%            subplot(2,2,3); plot(serieB); title('serieB'); grid on;
%            subplot(2,2,[2 4]); surface(g);ylabel('serieB'); xlabel('serieA'); axis image  ;
%            hold on; surf(pathArea); hold off;