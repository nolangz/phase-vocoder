%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         code for part two
%%%         Author:Nuolin Lai
%%%         Create Date:11/11/2020
%%%         Last modify date:11/25/2020
%%%         demo song 'mozart.wav'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear work space
clc;
clear all;

%import audio data
[x,Fs]  = audioread('mozart.wav');

%transfer to stereo audio
if size(x,2)==2
    x=0.5*x(:,1)+0.5*x(:,2);
end
x=0.5*sum(x,2);

%define frame time
t       = 100e-3;

%calculate frame size
N       = round(t*Fs);

%overlap factor
olf     =0.9;

%Analysis Hop size
HA      = round(N-N*olf);

%stretch factor
Q       = 1;

%if Q =1 Calculate N again
if Q==1
    N = round(HA/(1-olf));
end

%synthesis hop size
HS      = round(HA*Q);

%define NFFT
NFFT    =2^(floor(log2(N))+1);

%Qcreate hanning window
win     = 0.5*(1-cos(2*pi*(0:N-1)/N)).';

%calculate number of frame
L       = length(x);

%zeros padding in the front of x
x       = [zeros(N,1);x];

%Number of Frame
NF      = ceil((L+N)/HA);

%Total length of x after zeros padding in the front and end
L2      = (NF-1)*HA+N;

%Zeros padding in the end
x       = [x;zeros(L2-N-L,1)];

%create output vector y[n]
y       = zeros(HS*(NF-1)+N,1);

%create zeros vector
phi_m   = zeros(NFFT/2+1,1);
phi_m1  = zeros(NFFT/2+1,1);
omega_m1= zeros(NFFT/2+1,1);
theta_m = zeros(NFFT/2+1,1);

%frequency vector
wk      = (2*pi*(0:NFFT/2)/NFFT).';

%using for loop to create y[n]
for m=0:NF-1
    %window and FFT
    FFT         = fft(x(m*HA+1:m*HA+N).*win,NFFT);
    %Calculate Magnitude and phase of FFT 
    Xmag        = abs(FFT(1:NFFT/2+1));
    phi_m1      = angle(FFT(1:NFFT/2+1));
    %Calculate theta_m through emoga_m
    omega_m1    = wk+ppa(phi_m1-phi_m-wk*HA)/HA;
    theta_m     = theta_m+HS.*omega_m1;
    %reconstruct FFT
    Y           = Xmag.*exp(1i.*theta_m);
    %fix the error from pi by extracting the real part
    Y(1)        = real(Y(1));
    Y(NFFT/2+1) = real(Y(NFFT/2+1));
    % Hermitian-symmetric
    Y_m1        = [Y(1:NFFT/2+1);conj(Y(NFFT/2:-1:2))];
    %ifft
    frame_m1    = ifft(Y_m1,NFFT);
    %extracting the first N sample
    frame_m1    = frame_m1(1:N);
    %accumulation into y
    y(m*HS+1:m*HS+N,1) = frame_m1.*win+y(m*HS+1:m*HS+N,1);
    %Keep phi_m1 as phi_m to next loop
    phi_m       = phi_m1;
end

%gain factor 
gain = max(abs(x))/max(abs(y));

%reconstruct y
y=y*gain;

%play output audio
soundsc(y,Fs)

%create time vector for x
xt = ((0:length(x)-1)/Fs)';

%create time vector for x
yt = ((0:length(y)-1)/Fs)';
%plot x y and error and spectrogramwhen Q==1
if Q==1
    subplot(4,1,3);
    plot(xt,x,yt,y);
    subplot(4,1,4);
    plot(yt,x-y);
    subplot(4,1,1);
    MA1_S2119032_Lai_myspec(x,Fs,N,olf);
    subplot(4,1,2);
    MA1_S2119032_Lai_myspec(y,Fs,N,olf);
else
    %plot spectrogram
    subplot(2,1,1);
    MA1_S2119032_Lai_myspec(x,Fs,N,olf);
    subplot(2,1,2);
    MA1_S2119032_Lai_myspec(y,Fs,N,olf);
end