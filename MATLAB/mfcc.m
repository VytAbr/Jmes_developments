%%%% Mel-frequency cepstrum
% 1    Take the Fourier transform of (a windowed excerpt of) a signal.
% 2    Map the powers of the spectrum obtained above onto the mel scale, using triangular overlapping windows.
% 3    Take the logs of the powers at each of the mel frequencies.
% 4    Take the discrete cosine transform of the list of mel log powers, as if it were a signal.
% 5    The MFCCs are the amplitudes of the resulting spectrum.
%%%%
function [c FilterMask] = mfcc(signalas,c128) % signalas = 256 samples
% eile = 20;
NumOfFilter = 20;
dct_len     = 20;
FilterMask = mfcc_filt(NumOfFilter,50,5513,c128); % (order,lower_Freq,higher_Freq,len)

power = abs(fft(signalas)).^2;
power_mel = zeros(1,dct_len);
for i=1:dct_len
    power_mel(i) = sum(FilterMask(:,i).*power(1:c128)');
end
b = ones(1,dct_len);
for i=1:dct_len
    if power_mel(i) ~= 0
        b(i) = log2(power_mel(i));
    end
end

c = dct(fix(b),dct_len);
c = c(2:13); % grazina 12 MFCC koeficientu

% figure(1); subplot(2,2,1); bar(c); grid on; axis([0 13 -10 60]); title('MFCC');
%            subplot(2,2,2); bar(b); grid on; axis([0 21 0 50]); title('Log Mel Power');
%            subplot(2,2,3); plot(signalas); grid on; axis([0 256 -128 128]); title('Signal');
%            subplot(2,2,4); bar(power_mel); grid on; axis([0 21 -1 4000000]); title('Mel Power');



