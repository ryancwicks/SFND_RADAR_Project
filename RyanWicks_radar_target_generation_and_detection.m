clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3.0e8;
%% User Defined Range and Velocity of target
x_0 = 110; % m
v_0 = -50; % m/s
 
%% FMCW Waveform Generation
B = c/2.0/1.0; % Bandwidth = c/2/range resolution
Tchirp = 5.5 * 2 * 200 / c; % chirp time is about 5.5 time the max round trip time
slope = B/Tchirp; % Units of Hz/s 

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
dt = Nd*Tchirp/(Nr*Nd-1);
Fs = 1/dt;

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    %For each time stamp update the Range of the Target for constant velocity.
    x = x_0 + v_0 * t(i);
    
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2.0 * pi * (fc * t(i) + slope * t(i)^2.0 / 2.0 ) );
    t_offset = t(i) - 2.0 * x / c;
    Rx (i)  = cos(2.0 * pi * (fc * t_offset + slope * t_offset^2.0 / 2.0 ) );
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);
    
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix_2d = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
single_range = abs(fft(Mix_2d(:, 1))/Nr);
single_range = single_range(1:Nr/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output
f = Fs*(0:(Nr/2))/Nr;
r = f * c * Tchirp / 2 / B;
plot(r,single_range) 
title('Single-Sided Amplitude Spectrum of first FFT, plotted as range.')
xlabel('Range (m)')
ylabel('|P1(f)|')

axis([0, 200, 0, 0.3]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 4;
Gr = 5;
Gd = 2;

Ntot = (2*Tr + 2*Gr + 1)*(2*Td + 2*Gd + 1);
Nt = Ntot - (2*Gr + 1)*(2*Gd + 1);

% offset the threshold by SNR value in dB
offset = 7.0; %db

%Rather than looping, I set up a kernel based on the cell properties that
%can be convolved wth the signal to produce a noise array.
kernel = zeros(2*Tr + 2*Gr + 1, 2*Td + 2*Gd + 1) + 1.0/Nt;
kernel((Tr + 1):(Tr+2*Gr+1), (Td+1):(Td+2*Gd+1)) = 0.0;
pow_rdm = db2pow(RDM);

noise_level = conv2 (pow_rdm, kernel, "same");
noise_level = pow2db (noise_level) + offset;

RDM = int8((RDM-noise_level) > 0);
RDM(1:(Tr+Gr+1), :) = 0;
RDM((Nr/2 - Tr - Gr):Nr/2, :) = 0;
RDM(:, 1:(Td+Gd+1)) = 0;
RDM(:, (Nd - Td - Gd):Nd) = 0;

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM);
colorbar;


 
 