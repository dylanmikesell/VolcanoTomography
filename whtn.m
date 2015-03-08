% needs convsm1Df
function d1nw=whtn(flo,fhi,dt,d1n)


% whiten, then hanning taper back to 0.1-0.9 Hz?
% make a suitable hanning taper %
tlen = 361;
% tpr = hanning(tlen)';
% d1nf = d1nf.*[tpr(1:(tlen/2)) ones(1,(length(d1nf)-tlen)) tpr(((tlen/2)+1):tlen) ];
% d2nf = d2nf.*[tpr(1:(tlen/2)) ones(1,(length(d2nf)-tlen)) tpr(((tlen/2)+1):tlen) ];
tpr = hanning(tlen)';
%tpr2 = [ zeros(1,180) tpr(1:181) ones(1,((.9-.1)*36000*.1-1)) tpr(181:-1:1) zeros(1,(18000-180-181-2879-181)) ];
%tpr2 = [ zeros(1,180) tpr(1:181) ones(1,((.9-.1)*36000*.1-1)) tpr(181:-1:1) zeros(1,(18000-180-181-((.9-.1)*36000*.1-1)-181)) ];
tpr2 = [ zeros(1,180) tpr(1:181) ones(1,((fhi-flo)*36000*.1-1)) tpr(181:-1:1) zeros(1,(18000-180-181-((fhi-flo)*36000*.1-1)-181)) ];
tpr3 = [ 0 tpr2(18000:-1:2) tpr2 ];
% smooth in the fftshift mode and then hanning taper
% length of smoother is another knob to turn
% 180 means there is a window of .05 Hz on one side of arm
% i have to watch out for division by zero
% 9*10 is a one side window of .025 Hz
% smoothing
smlen = 90;
d1nff = fftshift(fft(d1n));
%size(d1nff./(convsm1Df(abs(d1nff)',smlen)+(eps*mean(abs(d1nff)))))
%size(tpr3)
d1nw = real(ifft(ifftshift((d1nff./(convsm1Df(abs(d1nff)',smlen)+(eps*mean(abs(d1nff))))).*tpr3)));
% no smoothing
% d1nff = fftshift(fft(d1n));
% d1nw = real(ifft(ifftshift((d1nff./(abs(d1nff)+(eps*mean(abs(d1nff))))).*tpr3)));
% d2nff = fftshift(fft(d2n));
% d2nw = real(ifft(ifftshift((d2nff./(abs(d2nff)+(eps*mean(abs(d2nff))))).*tpr3)));