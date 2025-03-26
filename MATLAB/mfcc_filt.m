function M = mfcc_filt(NumOfFilter,Flo,Fhi,len128)

% ############################## filter banks #############################
% NumOfFilter = 12;
Fhz0 = [Flo Fhi];
% Fhz0 = [50 5513]; % [300 8000]; % freq range
Fme0 = 1127*log(1+Fhz0/700); % Fme0 = 2595*log10(1+Fhz0/700);
dFme0 = Fme0(2) - Fme0(1);
dFme = dFme0/(NumOfFilter+1);
Fme = Fme0(1); % 1-st
for i=1:NumOfFilter+1
    Fme(i+1) = Fme0(1) + dFme*i; % 2-nd .. last
end
Fme = Fme';
Fhz = 700*(exp(Fme/1127)-1); % Fhz = 700*(10.^(Fme/2595)-1);
% ########################## x axis normalization #########################
norm_const = Fhz0(2)/len128;
Fhz_norm = round(Fhz/norm_const);
% ############################ banks formation ############################
BankMatrix = zeros(128,NumOfFilter);
for i=2:NumOfFilter+1
    for j=1:2
        if j==1
            d=Fhz_norm(i)-Fhz_norm(i-1);
            for k=1:d
                BankMatrix(Fhz_norm(i-1)+k,i-1) = k/d;
            end
        else
            d=Fhz_norm(i+1)-Fhz_norm(i);
            for k=1:d-1
                BankMatrix(Fhz_norm(i+1)-k,i-1) = k/d;
            end
        end
    end
end
M = BankMatrix;