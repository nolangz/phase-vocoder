%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         code for part three
%%%         Author:Nuolin Lai
%%%         Create Date:11/11/2020
%%%         Last modify date:11/25/2020
%%%         Define a function for Spectrogram plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define fuction
function myspec = MA1_S2119032_Lai_myspec(x,Fs,N,olf)
%time
t    = N/Fs;
%Analysis hop size
HA   = round(N-N*olf);
%Length of x
L    = length(x);
%Zeros padding in the front of x
x    = [zeros(N,1);x];
%number of frame
NF   = floor((L+N)/HA);
%length after zeros padding in the front and end
L2   = (NF-1)*HA+N;
%zeros padding in the end
x    = [x;zeros(L2-N-L,1)];
%hanning window
win  = 0.5*(1-cos(2*pi*(0:N-1)/N)).';
%number of NFFT
NFFT = 2^(floor(log2(N))+1);
%zeros matrix of STFT
STFT = zeros(NFFT,NF);
%for loop to fill the STFT matrix
for m=0:NF-1
    %accumulation
    STFT(:,m+1) = fft(x(m*HA+1:m*HA+N).*win,NFFT);
end
%take the first NFFT/2 number
STFT     = STFT(1:NFFT/2,:);
%time vector
x_time   = L/Fs.*(0:NF-1)/(NF-1);
%frequency vector
y_freq   = Fs/2*(0:NFFT/2-1)/(NFFT/2-1);
%magnitude of STFT
STFT_mag = 20*log10(abs(STFT));
%plot the image
image(x_time,y_freq,STFT_mag,'CDataMapping','scaled');
%set y scale to start 0 from bottom
set(gca,'YDir','normal');
%create colorbar
colorbar;
%set colormap limits
lim      = caxis;
caxis([max(STFT_mag(:))-60 max(STFT_mag(:))]);
%xlabel
xlabel('Time (s)');
%ylabel
ylabel('Frequency (Hz)');
title({'Spectrogram with',['frame Length=',num2str(N/Fs),' second'],['over lap factor=',num2str(olf*100),'%']});
end