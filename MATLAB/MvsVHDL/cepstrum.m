function C = cepstrum(s,eile)

s = s.*hann(256)';
a = abs(fft(s));
b = ones(1,256);
for i=1:256
    if a(i) ~= 0
        b(i) = a(i);
    end
end
c = fix(real(fft(fix(log2(b)))));
% c = cceps(s);
C = abs(c(2:eile+1));