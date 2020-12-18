clear all
clc
figure(1);clf
%y=[ 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1]; % spike train
SPIKE_TRAIN_LENGTH = 100;
y=double(logical(rand(SPIKE_TRAIN_LENGTH,1)<0.5));
plot(y)
hold on
sigma=1;
x=[-100:sigma:100];

%k = exp(-(x/sigma).^2/2)/(sigma*sqrt(2*pi)); % Gaussian kernel.
k = gausswin(3);
z=conv(y,k);
z=z(ceil(length(k)/2):end-floor(length(k)/2));
h=plot(z, '-r');
h.Color(4) = 0.5;  % 50% transparent

legend('Spike Train', 'Rate');