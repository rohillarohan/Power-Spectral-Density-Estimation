%% Power Spectral Density Estimation
clc
clearvars
close all
%% Generating Colored Noise
fs = 10e3;
ts = 1/fs;
t1 = 0:ts:500*ts-ts;
h1 = tf([1 0 0],[1 -1.585 0.96],1);
whitNos = randn([1,500]);
disp('Varaince of the White Noise')
disp(var(whitNos))
[h1,w1] = impz([1 0 0],[1 -1.585 0.96]);
[H1,W1] = freqz(h1,fs);

figure;
subplot(221)
plot((0:length(whitNos)-1),whitNos,'linewidth',1)
title('white Noise')
xlabel('Time Index')
ylabel('Value')
axis([0 length(whitNos)-1 -4 3])
grid

subplot(222)
plot(w1,h1,'linewidth',1)
title('Immpulse Response of the Filter')
xlabel('Time Index')
ylabel('Value')
axis([0 length(h1)-1 -1.5 2])
grid

subplot(223)
plot(W1,20*log10(abs(H1)),'linewidth',1)
title('Magnitude Response of the Colouring Filter (Ideal Spectrum)')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
axis([0 max(W1) -95 -45])
grid

clrdNos = conv(whitNos,h1);
[CN,WCN] = freqz(clrdNos,fs);
subplot(224)
plot(WCN,20*log10(abs(CN)),'linewidth',1)
title('Magnitude Response of the Coloured Noise')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
axis([0 max(WCN) -110 -10])
grid

%% Periodogram with 1 Realization
[p1ClrdNos,pW1] = periodogram(clrdNos);

figure;
plot(W1,20*log10(abs(H1)),'g','linewidth',1)
hold on 
plot(pW1,20*log10(p1ClrdNos),'r','linewidth',1)
grid
title('Response of Regular Periodogram vs Ideal Response')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
legend('Ideal Response','Response of Regular Periodogram')
axis([0 pi -100 60])

disp('Variance of 1 Realization using regular Periodogram')
disp(var(p1ClrdNos))
%% Periodogram with 50 different realizations
nwWhitNos = zeros(50,500);

for n = 1:50
    nwWhitNos(n,:) = randn([1,500]);
    nwClrdNos(n,:) = conv(nwWhitNos(n,:),h1);
    [p2ClrdNos(n,:),pW2] = periodogram(nwClrdNos(n,:));
end
ensmblAvg = mean(p2ClrdNos);
pClrdNosdB=20*log10(p2ClrdNos');

figure;
subplot(211)
plot(pW2,pClrdNosdB)
grid
title('Periodogram Display of 50 Different Realizations')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
axis([0 pi -150 70])

ensmblAvgdB=20*log10(ensmblAvg);
subplot(212)
plot(pW2,ensmblAvgdB,'b','linewidth',1)
hold on
plot(W1,20*log10(abs(H1)),'g','linewidth',1)
hold on 
plot(pW2,20*log10(p2ClrdNos(n,:)),'r','linewidth',1)
grid
title('Average of 50 Different Realizations vs Ideal Response')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
legend('Ensemble Average','Ideal Response','Response of 1 Realization')
axis([0 pi -100 60])

disp('Variance of 1 Realization using regular Periodogram')
disp(var(p2ClrdNos(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations using regular Periodogram')
disp(var(ensmblAvg))

% Variance of the Ensemble Average vs 1 Realization:
% We see that the variance of one realization is greater than that of the
% ensemble average and is really obvious. We are averaging different
% variables, which leads to a lower variance. The local bias of the single
% frequency spectrum is significantly lower and the high frequency bias is
% significantly higher.

%% Blackman Tukey Method
wind = hamming(984-1);

for n = 1:50
    [btClrdNosPer(n,:),BTW]=per_smooth(nwClrdNos(n,:),wind,length(nwClrdNos(n,:))/2);
end
btClrdNosPerAve=mean(btClrdNosPer);
figure;
btClrdNosPerdB=20*log10(btClrdNosPer');
subplot(2,1,1)
plot(BTW,btClrdNosPerdB)
grid
title('Blackma-Tukey with Hamming Window')
xlabel('Normalized Frequency')
ylabel('LogMagnitude');
axis([0 pi -80 10+max(max(btClrdNosPerdB))])

subplot(2,1,2)
btClrdNosPerAvedB=20*log10(btClrdNosPerAve);
plot(BTW,btClrdNosPerAvedB,'b','linewidth',1)
hold on
plot(W1,20*log10(abs(H1)),'g','linewidth',1)
hold on
plot(BTW,20*log10(btClrdNosPer(n,:)),'r','linewidth',1)
grid on
axis([0 pi -100 80])
title('Average of 50 Different Realizations')
xlabel('Normalized Frequency')
ylabel('Log Magnitude')
legend('Average of 50 Different Realizations','Ideal Response','Response of 1 Realization')

disp('Variance of 1 Realization Using Blackman-Tukey')
disp(var(btClrdNosPer(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Blackman-Tukey')
disp(var(btClrdNosPerAve))

% Blackman-Tukey Method:
% What Balckman-Tukey does is just takes the periodogram and smooths it using
% a window. This window could be any window, we have used Hamming window
% here. After applying the window the abberations in the signle realization
% spectrum have reduced. Similarly the ensemble average spectrum is also
% very smooth. 
%% MVDR with 200 Tap FIR Filter
w = 0:pi/1024:pi;
for n = 1:50
    [mv1ClrdNosSpctr(n,:)]=MV(nwClrdNos(n,:),199,w);
end
figure;
subplot(211)
mv1ClrdNosSpctrdB = 10*log10(mv1ClrdNosSpctr');
plot(w/(2*pi),mv1ClrdNosSpctrdB)
grid
axis([0 1/2 min(min(mv1ClrdNosSpctrdB))-10 10+max(max(mv1ClrdNosSpctrdB))])
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Minimum Variance Method (200 Length Filters)')

subplot(212)
plot(w/(2*pi),10*log10(mv1ClrdNosSpctr(n,:)'),'r','linewidth',1)
grid
hold on
mv1ClrdNosSpctrAve=mean(mv1ClrdNosSpctr);
mv1ClrdNosSpctrAvedB=20*log10(mv1ClrdNosSpctrAve);
plot(w/(2*pi),mv1ClrdNosSpctrAvedB,'b','linewidth',1)
hold on
plot(W1/(2*pi),20*log10(abs(H1)),'g','linewidth',1)
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Minimum Variance Spectrum Comparision (200 Length Filters)')
legend('1 Realization','Ensemble Average','Ideal Response')

disp('Variance of 1 Realization Using Minimum Variance (200 Length Filters)')
disp(var(mv1ClrdNosSpctr(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Minimum Variance (200 Length Filters)')
disp(var(mv1ClrdNosSpctrAve))

% Minimum Variance Method (Length 200 FIR Filter):
% Just by looking at the plot one can say that 200th order local FIR
% filters aren't enough. The one realization spectrum is not even close to
% the ensemble average spectrum. In addition to that the variation are so
% huge that it is difficult to find the main frequency. The ensemble
% average spectrum is not good either. for sure it is smooth but lower high
% frequency bias gives rise to higher local bias, and that is quite visible
% here.
%% MVDR with 500 tap FIR Filters 

for n = 1:50
    [mv2ClrdNosSpctr(n,:)]=MV(nwClrdNos(n,:),499,w);
end
figure;
subplot(211)
mv2ClrdNosSpctrdB = 10*log10(mv2ClrdNosSpctr');
plot(w/(2*pi),mv2ClrdNosSpctrdB)
grid
axis([0 1/2 min(min(mv2ClrdNosSpctrdB))-10 10+max(max(mv2ClrdNosSpctrdB))])
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Minimum Variance Method (500 Length Filters)')

subplot(212)
plot(w/(2*pi),20*log10(mv2ClrdNosSpctr(n,:)'),'r','linewidth',1)
grid
hold on
mv2ClrdNosSpctrAve=mean(mv2ClrdNosSpctr);
mv2ClrdNosSpctrAvedB=20*log10(mv2ClrdNosSpctrAve);
plot(w/(2*pi),mv2ClrdNosSpctrAvedB,'b','linewidth',1)
hold on
plot(W1/(2*pi),20*log10(abs(H1)),'g','linewidth',1)
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Minimum Variance Spectrum Comparision (500 Length Filters)')
legend('1 Realization','Ensemble Average','Ideal Response')

disp('Variance of 1 Realization Using Minimum Variance (500 Length Filters)')
disp(var(mv2ClrdNosSpctr(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Minimum Variance (500 Length Filters)')
disp(var(mv2ClrdNosSpctrAve))

% Minimum Variance Method (Length 500 FIR filters):
% Now the difference can be noted even by a non signal processing
% engineer!! 500 length filters have done the job quite well. The single
% realization spectrum matches the shape of the ensemble average spectrum
% (Not completely though but still has some improvement). The main
% frequency can be spotted from the spectrum. Even the variation through
% out the spectrum have lower magnitude but now occur at a higher rate. The
% ensemble average spectrum is sharper as compared to the previous attempt.
% But as said earlier, lower local bias comes at a cost of higher high
% frequency bias.
%% Welch's Method with no Segment Overlaps
N = 984; % Number of samples in the signal
K = 3; % Number of segments
L1 = N/K; % In case ofNo overlap
wind = hamming(L1); % Rectangular Window

% Generating 50 sample signals and their periodograms
for n=1:50
    [wlch1ClrdNosSpctr(n,:),Wlch]=pwelch(nwClrdNos(n,:),wind,L1/2,2^12);
end
wlch1ClrdNosSpctrAve=mean(wlch1ClrdNosSpctr);

figure;
wlch1ClrdNosSpctrdB=20*log10(wlch1ClrdNosSpctr');
subplot(2,1,1)
plot(Wlch/(2*pi),wlch1ClrdNosSpctrdB)
axis([0 1/2 -80 10+max(max(wlch1ClrdNosSpctrdB))])
title('Welch Method (No Overlap)')
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
grid

wlch1ClrdNosSpctrAvedB=20*log10(wlch1ClrdNosSpctrAve);
subplot(2,1,2)
plot(Wlch/(2*pi),20*log10(wlch1ClrdNosSpctr(n,:)'),'r','linewidth',1)
hold on
plot(Wlch/(2*pi),wlch1ClrdNosSpctrAvedB,'b','linewidth',1)
hold on
plot(W1/(2*pi),20*log10(abs(H1)),'g','linewidth',1)
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Welch Method Spectrum Comparision (No Overlap)')
legend('1 Realization','Ensemble Average','Ideal Response')
axis([0 1/2 -100 10+max(wlch1ClrdNosSpctrAvedB)])
grid

disp('Variance of 1 Realization Using Welch Method (No Overlap)')
disp(var(wlch1ClrdNosSpctr(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Welch Method (No Overlap)')
disp(var(wlch1ClrdNosSpctrAve))

% Welch's Method with No Overlap:
% Using the Welch's method with no overlap the spectrum of just one
% realization has abberations like ususal but a point to note is that the
% local bias around the main frequency is high. The ensemble average
% spectrum is comparitively smooth.
%% Welch's Method with 50% Overlap Between Adjacent Segments

N=984; % Number of samples in the signal
K=3; % Number of segments
L2=(2*N)/(K+1); % In case of 50% overlap
wind=hamming(L2); % Rectangular Window

for n=1:50
    [wlch2ClrdNosSpctr(n,:),Wlch]=pwelch(nwClrdNos(n,:),wind,L2/2,2^12);
end
wlch2ClrdNosSpctrAve=mean(wlch2ClrdNosSpctr);

figure;
wlch2ClrdNosSpctrdB=20*log10(wlch2ClrdNosSpctr');
subplot(2,1,1)
plot(Wlch/(2*pi),wlch2ClrdNosSpctrdB)
axis([0 1/2 -80 10+max(max(wlch2ClrdNosSpctrdB))])
title('Welch Method (50% Overlap)')
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
grid

wlch2ClrdNosSpctrAvedB=20*log10(wlch2ClrdNosSpctrAve);
subplot(2,1,2)
plot(Wlch/(2*pi),20*log10(wlch2ClrdNosSpctr(n,:)'),'r','linewidth',1)
hold on
plot(Wlch/(2*pi),wlch2ClrdNosSpctrAvedB,'b','linewidth',1)
hold on
plot(W1/(2*pi),20*log10(abs(H1)),'g','linewidth',1)
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Welch Method Spectrum Comparision (50% Overlap)')
legend('1 Realization','Ensemble Average','Ideal Response')
axis([0 1/2 -100 10+max(wlch2ClrdNosSpctrAvedB)])
grid

disp('Variance of 1 Realization Using Welch Method (50% Overlap)')
disp(var(wlch2ClrdNosSpctr(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Welch Method (50% Overlap)')
disp(var(wlch2ClrdNosSpctrAve))

% Welch's Method with 50% overlap:
% On comparing the response with the one with no overlap we can see that
% the spectrum of just one realization has lower magnitude of variations
% but this variations are now at a quicker rate. The ensemble average
% spectrum here has a lower local bias and higher (comparitively) high frequency
% bias.
%% Multitapere Method
N = 984;
dpss_params = {2.5,3};
for n=1:50
    [mtpr1ClrdNosSpctr(n,:) MT] = pmtm(nwClrdNos(n,:),dpss_params);
end
mtpr1ClrdNosSpctrAve=mean(mtpr1ClrdNosSpctr);

figure;
mtpr1ClrdNosSpctrdB=20*log10(mtpr1ClrdNosSpctr');
subplot(2,1,1)
plot(MT/(2*pi),mtpr1ClrdNosSpctrdB)
axis([0 1/2 -80 10+max(max(mtpr1ClrdNosSpctrdB))])
title('Multitapered Method (DPSS window)')
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
grid

mtpr1ClrdNosSpctrAvedB=20*log10(mtpr1ClrdNosSpctrAve);
subplot(2,1,2)
plot(MT/(2*pi),20*log10(mtpr1ClrdNosSpctr(n,:)'),'r','linewidth',1)
hold on
plot(MT/(2*pi),mtpr1ClrdNosSpctrAvedB,'b','linewidth',1)
hold on
plot(W1/(2*pi),20*log10(abs(H1)),'g','linewidth',1)
xlabel('Normalized Frequency')
ylabel('Log Magnitude');
title('Multi Taper Spectrum Comparision (DPSS Window)')
legend('1 Realization','Ensemble Average','Ideal Response')
axis([0 1/2 -100 10+max(mtpr1ClrdNosSpctrAvedB)])
grid

disp('Variance of 1 Realization Using Multi Taper Method with DPSS window')
disp(var(mtpr1ClrdNosSpctr(n,:)))
disp('Variance of the Ensemble Average of 50 Realizations Using Multi Taper Method with DPSS window')
disp(var(mtpr1ClrdNosSpctrAve))

% Multitaper:
% From the plots we can see that the spectrum of the ensemble average is really sharp at the main
% frequency and has abberations through out the band. In case of one
% realization it seems to have very high variations as compared to the
% ensemble average plot. Talking about Bais, One single realization seems
% to have a slightly narrow spectrum at the main frequency. It has Lower
% Local bias and higer high frequency bias
%% Discussion and Comparision of the Algorithms
% Different methods gave us different spectrum. Some had smoother spectrum
% like Multitaper/Welch (No Overlap), whereas some had sharper spectrum with 
% high resolution like Blackman-Tukey. We also saw that some methods gave us 
% comparitively good responses even with a single realizations like Multitaper 
% when compared to Minimum Variance (200/500 length). Now depending on the 
% application and the amount of money a company can spend to make a sophisticated 
% product decides on to which method will be pciked. If A company can have 
% different realizations they can definitely go for Welch's Method (as it is really 
% smooth), whereas it is highly likely that only one realization will be available. 
% In that case one will go for Multitaperd as it gives a workable spectrum only 
% with one realization. Sometimes what matters is resolution to seperate
% nearby frequencies. In this case Minimum Variance does the trick. It can
% be seen that it has a really nice peak. And if we have a good band pass
% filter then this method would be enough. SImilarly if resolution is
% needed and only one realization is available one can even go for a
% periodogram but it would need a really good narrow band filter to take
% the peak out.