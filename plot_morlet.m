%function [ output_args ] = plot_windowed_fft( input_args )

%%
T = 3600;
step = 100;
k=floor(log(T)/log(1.5));
% widths = floor(1.5.^[1:k]);
widths = 1:step:(2*T);  %% 2T!!
L = length(data);
mat = abs(mycwt(data, widths));

f=widths/T;

figure
% colormap(hot)
% imagesc(time(Ts)./T, f(1:flim), mat(1:flim, :))
imagesc(time./T, f, mat)
title('Morlet Transform')
xlabel('Revolution')
ylabel ('period (revolution)') %ylabel('frequency') 
set(gca,'Ydir','normal')

set(gca, 'XTick', [10:29])

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
% figure 
% norm = sum(sqrt( mat(:,2:end)-mat(:,1:(end-1)).^2 ));  % L2 norm
% plot(time(Ts)./T, norm)