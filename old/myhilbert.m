function myhilbert()

n = 21;
rnums = randn(n,1);

f = fft(rnums);

complexf = 1i*f;

posF = 2:floor(n/2)+mod(n,2);
negF = ceil(n/2)+1+~mod(n,2):n;

f(posF) = f(posF) + -1i*complexf(posF);
f(negF) = f(negF) + 1i*complexf(negF);

hilbertx = ifft(f);
hilbertm = hilbert(rnums);

figure; plot(abs(hilbertx)); hold on; plot(abs(hilbertx),'o')
figure; plot(angle(hilbertx)); hold on; plot(angle(hilbertm),'o')

end 