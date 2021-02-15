%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         code for part one
%%%         Author:Nuolin Lai
%%%         Create Date:11/11/2020
%%%         Last modify date:11/25/2020
%%%         demo song 'mozart.wav'
%%%         when Q = 1.25 N=4410 HA=1764 is a suitable choice, 
%%%         overlap factor= 0.55
%%%         when Q = 0.75 N=4410 HA=1985 is a suitable choice, 
%%%         overlap factor=0.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear work space
clc;
clear all;

%import audio data
[x,Fs] = audioread('mozart.wav');

%transfer to stereo audio
if size(x,2)==2
    x=0.5*x(:,1)+0.5*x(:,2);
end
x = 0.5*sum(x,2);

%define frame time
t      = 100e-3;

%calculate frame size
N      = round(t*Fs);

%overlap factor
olf    = 0.75;

%Analysis Hop size
HA     = round(N-N*olf);

%stretch factor
Q      = 1;

%if Q =1 Calculate N again
if Q==1
    N = HA/(1-olf);
end

%Synthesis hop size
HS     = round(HA*Q);

%Qcreate hanning window
win    = 0.5*(1-cos(2*pi*(0:N-1)/N)).';

%calculate number of frame
L      = length(x);

%zeros padding in the front of x
x      = [zeros(N,1);x];

%Number of Frame
NF     = ceil((L+N)/HA);

%total length
L2     = (NF-1)*HA+N;

%Zeros padding in the end
x      = [x;zeros(L2-N-L,1)];

%create zeros output vector y[n]
y      = zeros(HS*(NF-1)+N,1);

%using for loop to create y[n]
for m=1:NF
    %accumulation into y
    y((m-1)*HS+1:(m-1)*HS+N,1) = x((m-1)*HA+1:(m-1)*HA+N,1).*win+y((m-1)*HS+1:(m-1)*HS+N,1);
end

%play audio
soundsc(y,Fs);

%gain factor 
gain = max(abs(x))/max(abs(y));

%reconstruct y
y=y*gain;

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