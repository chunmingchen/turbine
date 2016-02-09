%function [ output_args ] = plot_windowed_fft( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% Fs = L;            % Sampling frequency
% f = Fs*(0:(L/2))/L/(L/3600);
% % t = 1./f(2:end);

%%
FROM_PROBES = 1
if FROM_PROBES
    T = 3600
    W = T % window
    step = 10 % 10
else
    T = 144
    W =  T
    step = 1
end
L = length(data)
Ts = step:step:L-W;
mat = zeros(floor(W/2)+1, length(Ts) );
matp = mat; % phase
mat1 = zeros(floor(W/2)+1, length(Ts) );
count = 1;
for t=Ts
    %display(sprintf('%d %d', ((t-1)*W+1) ,(t*W)))
    Y = fft(data(1, t:(t+W-1) ))';
    P2 = abs(Y)/W; % Fourier amplitude
    P1 = P2(1:W/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    if FROM_PROBES
        P1 = min(P1, 0.01); % clamp !!!!! Remove this for any new signal
    end
    
    count = count +1;
    mat(:,count) = P1;
    
%     P3 = angle(Y)/pi*180; % Fourier phase (degree)
%     P3(find(P2<0.0001))=0;
%     P3 = P3(1:W/2+1);
%     matp(:, count) = P3;
    
%     
%     plot(P1)
%     pause
end
%f = (0:20:(L/2))/(L/T);
%f = (0:4:(L/2*W/T))/(L/T*W/T);
f = (0:(W/2))/(W/T);

if FROM_PROBES
    flim = 100;
else
    flim = size(mat, 1);
end

figure
% colormap(hot)
imagesc(time(Ts)./T, f(1:flim), mat(1:flim, :))
title('windowed fft')
xlabel('Revolution')
ylabel ('Frequency (1/revolution)') %ylabel('frequency') 
set(gca,'Ydir','normal')

% set(gca, 'XTick', [10:29])

%surface(1:size(mat,1), 1:L-W,  mat)

% colormap(hot)
% cmap = colormap(hot);
% cmap = cmap([1:44,64],:);
% cmap = flipud(cmap);
% colormap(cmap)
colormap parula
colorbar

% figure
% imagesc(time(Ts)./T, f(1:flim), matp(1:flim, :))
% title('windowed fft phase')
% xlabel('Revolution')
% ylabel('Phase (degree)')
% set(gca, 'Ydir', 'normal')

% show difference
figure 
diff =  mat(:,2:end)-mat(:,1:(end-1));
norm = sqrt(sum(diff.^2 , 1)) ;  % L2 norm
plot(time(Ts)./T, norm)
