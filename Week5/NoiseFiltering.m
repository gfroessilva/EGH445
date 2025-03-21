%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Digital Control
% **** Main Functions
%     **** De-noising Filter Design Resources
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------



%% Filter (Anti-Alias) Resources  (requires update)
% Noise Filter (Low Pass FIR)
% Design a minimum-order lowpass FIR filter with normalized passband 
% frequency 0.25pi rad/s, stopband frequency 0.35pi rad/s, passband ripple 0.5 dB, and 
% stopband attenuation 65 dB. Use a Kaiser window to design the filter. 
lpFilt = designfilt('lowpassfir','PassbandFrequency',0.25, ...
         'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'DesignMethod','kaiserwin');
fvtool(lpFilt)                  % Visualize Filter Magnitude Response.  
                            
dataIn = rand(1000,1);          % Vector Random Data
dataOut = filter(lpFilt,dataIn);% Filter Random Data.
figure(9); hold on
plot(dataIn,'r')
plot(dataOut,'b')


% Nose Filter (Low Pass IIR)
% Design a lowpass IIR filter with order 8, passband frequency 150 Hz, and 
% passband ripple 0.2 dB. Specify a sample rate of 1000 Hz. Visualize the 
% magnitude response of the filter. Use it to filter a 1000-sample random signal.
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',150,'PassbandRipple',0.5, ...
         'SampleRate',450);
fvtool(lpFilt)
dataIn = randn(1000,1);
dataOut = filter(lpFilt,dataIn);

% Frequncy Response from Transfer Function
b0 = 0.05634;
b1 = [1  1];
b2 = [1 -1.0166 1];
a1 = [1 -0.683];
a2 = [1 -1.4461 0.7957];
b = b0*conv(b1,b2);
a = conv(a1,a2);

[h,w] = freqz(b,a,'whole',2001);
plot(w/pi,20*log10(abs(h)))
ax = gca;
ax.YLim = [-100 20];
ax.XTick = 0:.5:2;
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% Help 
% https://au.mathworks.com/help/signal/ref/designfilt.html
% https://au.mathworks.com/help/signal/ref/digitalfilter-class.html