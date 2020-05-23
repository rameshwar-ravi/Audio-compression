%Matlab code for signal compression using DCT
clear all
close all
clc
%%input the audio
%ff-16b-1c-44100hz.wav
[a, fs] = audioread('choose.wav');
prev = 1;
data = 0;
a = a(1:160000);
subplot(2,1,1)
plot(a)
title("Original Signal")
%input audio

sound(a, fs);
pause(5)


%Since we've total 40,000 samples and we assume they all represented by
%different symbols, so total symbols will be 2^n, where n = log2(length(a))
power = ceil(log2(length(a)))
input_signal_bits = length(a)*power
c = 0;


ppo = 0;
for i = 255:255:length(a)
    k = dct(a(prev:i));
    prev=prev+255;
    ppo = ppo+1;
    data = [data transpose(k(1:128))];
end

%plot(data)
%title("Compressed DCT of Signal")

%% Quantization
%now do quantization
%partition the data
%get the max and min
data_max = max(data)
data_min = min(data)
%define the step size 
step_size = 0.001;
%partition into the levels
levels = data_min:step_size:data_max;
%define the values 
values = data_min-(step_size)/2:step_size:data_max+(step_size)/2;
[index, quants] = quantiz(data, levels, values); %quants is my quantized data

%% Huffman encoding
%computing the probability distribution
freq = zeros(1,length(values));
for i = 1:length(quants)
    freq(round((quants(i)-(data_min-step_size/2))*100) + 1) = freq(round((quants(i)-(data_min-step_size/2))*100) +1) + 1;
end
%convert it into prob
prob = freq/length(quants);
%huffman encoding
[dict, avglen] = huffmandict(1:length(prob),prob);
%now count the number of bits 
tbits = 0;
for i = 1:length(quants)
    %round((quants(i)-(data_min-step_size/2))*100) + 1 => will give me the
    %dict index
    cell_i = dict(length(dict(round((quants(i)-(data_min-step_size/2))*100) + 1)),2);
    cell_mat = cell2mat(cell_i);
    tbits = tbits + length(cell_mat);
end
tbits

% now we have decompress the signal for our useage... 
% we have to do the idct and compute mse
prev = 1;
temp = [];

j = 0;



for i = 128:128:length(quants)
    
    k = [data(prev:i) zeros(1,127)];
    prev=prev+128;
    j =j + 1;
    temp  = [temp idct(k)]; 
end
diff =  length(a)-length(temp);
%padding with zeros to calc. MSE and PSNR
temp = [temp zeros(1, diff)];
psnr_ = psnr(a', temp)
mse_ = mse(a', temp)
snr_ = snr(a', temp)
sound(temp, fs)
pause(10)
subplot(2,1,2)
plot(temp)
title("Received Signal")

%Space saved
%calculate number of bits in original signal, after dct, after
%quantization, after huffman coding and bits of compressed signal
%compression_ratio_percentage = (tbits/input_signal_bits)*100
Compression_ratio = input_signal_bits/tbits

 


