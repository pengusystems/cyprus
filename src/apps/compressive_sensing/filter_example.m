clc;

% Parameters.
n_cs_samples = 50;
fs = 1000;
n_samples = fs + 1;
t = (0:1/fs:1)';
f = (-fs/2:fs/2)';
inputdata1 = sin(2*pi*20*t);
inputdata2 = sin(2*pi*130*t);

% Filter design.
% Fp - filter pass
% Fst - filter stop
% Ap - Allowed ripple around 0dB
% Ast - Attenuation in stop band
d = fdesign.lowpass('Fp,Fst,Ap,Ast',100,120,0.5,50,fs);
Hd = design(d);

% Filter and plot.
output1 = filtfilt(Hd.Numerator,1,inputdata1);
output2 = filtfilt(Hd.Numerator,1,inputdata2);
figure; plot(output1);
figure; plot(output2);