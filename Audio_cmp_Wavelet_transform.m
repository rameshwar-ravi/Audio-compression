clear;
close all;
clc;
filename = 'choose.wav';
[x,fs] =audioread(filename);
framesize=2048;
s=dir(filename);
input_size=s.bytes;
prcessed_out=0;
prcessed_out_length=zeros(7,1);
levels=5;
total_steps=ceil(length(x)/framesize);
decodedoutput=1;
codes_out=0;
entropy_length=0;
gain=0;
compressed_out=0;
storagegain=0;
total_len=0;
storage_len=0;
for i=1:1:total_steps
    %take care of remainder at the last dataframe
    bits=11;
    if(i~=total_steps)
        frame_data=x([framesize*(i-1)+1:framesize*i]);
    else
        frame_data=x([framesize*(i-1)+1:end]);
        bits=8;
    end
    %initial conditions for transform compression
    [c,l] = wavedec(frame_data,levels,'haar');
    %Wavelet decomposition vector, returned as a real-valued vector. 
    %The bookkeeping vector l contains the number of coefficients by level.
    [default1,default2,default3] = ddencmp('cmp','wv',frame_data);%default values
    [compressed,final_c,final_l,PERF0,PERFL2] = wdencmp('gbl',c,l,'haar',levels,default1,default2,default3);
    
    %compander Mu-law
    compand_out=compand(final_c,255,max(final_c),'mu/compressor');
    %quantizer
    delta=(max(compand_out)-min(compand_out))/2^8;
    intervals=[min(compand_out):delta:max(compand_out)];
    codes=[1 min(compand_out):delta:max(compand_out)];
    codes_out=[codes_out codes];
    [qindex,quantize_out]=quantiz(compand_out,intervals,codes);
    index=find(compand_out==0);
    adjust=quantize_out(index(1));
    quantize_out=quantize_out-quantize_out(index(1));
    prcessed_out=[prcessed_out quantize_out];
    %prcessed_out_length=[prcessed_out_length final_l];
    %extraxtion
    pdf = zeros(1,length(codes));
    for myvar = 1:1:length(codes)
        pdf(myvar)=(sum(quantize_out==codes(myvar)-adjust))/length(quantize_out);    
    end
    symbols=1:length(codes);
    [dict, average_len] = huffmandict(symbols,pdf);
   gain=gain+(length(frame_data)*bits)-(average_len*length(frame_data)); 
   storagegain=storagegain+(length(frame_data)*bits-average_len*258);
   total_len=total_len+length(frame_data)*bits;
   storage_len=storage_len+average_len*258;
    outframe=quantize_out;
    if(max(outframe)~=0)
        outframe=compand(outframe,255,max(outframe),'mu/expander');
    end
    out=waverec(outframe,final_l,'haar');
    decodedoutput=[decodedoutput out];
end
decodedoutput=decodedoutput(2:length(decodedoutput));
audiowrite('out.wav',decodedoutput,fs);
[xout,fs]=audioread('out.wav');
 sound(x,fs);
 pause(10);
sound(xout,fs);
time=0:1/fs:length(x)/fs-1/fs;
subplot(2,1,1);
plot(x);
subplot(2,1,2);
plot(xout,'r');
figure;
plot(x);title('Difference of Input Vs Output');
hold on
plot(xout,'r')
plot(x-xout,'y')
mse=mse(x,xout);
PSNR = psnr(x,xout);
se1=input_size/1024;
s11=dir('out.wav');
se2=s11.bytes/1024;
memory_in=184368;
memory_out=23479;
save=memory_out/memory_in;
disp('save in terms of sample');disp(save);
disp('gain in terms of bits due to huffman coding');disp(gain);
disp('gain in terms of bits due compression');disp(storagegain);
disp('input bits/output bit');
disp(total_len/storage_len);