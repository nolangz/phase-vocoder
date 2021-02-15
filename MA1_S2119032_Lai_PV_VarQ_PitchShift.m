%code for varible and pitch shift
%Author:Nuolin Lai
%Create Date:11/18/2020
%Last modify date:11/26/2020

%clear work space
clc;
clear all;

%import audio data
[x,Fs] = audioread('mozart.wav');

%cent
c       = -12;

%pitch shift ratio
ps      = 2^(c/12);

%transfer to stereo audio
if size(x,2)==2
    x=0.5*x(:,1)+0.5*x(:,2);
end
x=0.5*sum(x,2);

%define frame time
t = 20e-3;

%calculate frame size
N = round(t*Fs);

%overlap factor
olf=0.75;

%Analysis Hop size
HA= round(N-N*olf);

% input Q linear start
start=0.2;

% input Q linear end
stop=5;

%define number of FFT
NFFT=2^(floor(log2(N))+1);

%create hanning window
win = 0.5*(1-cos(2*pi*(0:N-1)/N)).';

%calculate number of frame
L = length(x);

%zeros padding in the front of x
x = [zeros(N,1);x];

%Number of Frame
NF = ceil((L+N)/HA);

%Total length of x after all processing
L2 = (NF-1)*HA+N;

%Zeros padding in the end
x = [x;zeros(L2-N-L,1)];

%Q vector
Q=linspace(start,stop,NF-1);

%HS vector
HS = round(HA.*Q);

%RS vector 
RS = round(HS.*ps);

%create zeros output vector y[n]
HS_sum = sum(RS);
Ly = HS_sum+N;
y = zeros(Ly,1);

%create zeros vector
phi_m     = zeros(NFFT/2+1,1);
phi_m1    = zeros(NFFT/2+1,1);
omega_m1  = zeros(NFFT/2+1,1);
theta_m   = zeros(NFFT/2+1,1);

%frequency vector
wk        = (2*pi*(0:NFFT/2)/NFFT).';

%using for loop to create y[n]
for m=1:NF-1
    %window and FFT
    FFT         = fft(x(m*HA+1:m*HA+N).*win,NFFT);
    %Calculate Magnitude and phase of FFT 
    Xmag        = abs(FFT(1:NFFT/2+1));
    phi_m1      = angle(FFT(1:NFFT/2+1));
    %Calculate theta_m through emoga_m
    omega_m1    = wk+ppa(phi_m1-phi_m-wk*HA)/HA;
    theta_m     = theta_m+RS(m).*omega_m1;
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
    y(sum(RS(1:m))+1:sum(RS(1:m))+N,1) = frame_m1.*win+y(sum(RS(1:m))+1:sum(RS(1:m))+N,1);
    %Keep phi_m1 as phi_m to next loop
    phi_m = phi_m1;
end
%resampling
y=interp1(y,(1:ps:Ly));

%play audio
soundsc(y,Fs);
