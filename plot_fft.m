idx = 5483896
i = mod(idx , nzcols)+1
t = floor(idx/nzcols)+1
selected = coeffs(:,i:nzcols:end);
selected = selected + repmat(freq_mean, 1, size(selected,2));
selected(1,:) = [];

figure;
subplot(3,3,1)
imagesc(selected);
subplot(3,3,2)
mesh(selected)
subplot(3,3,3)
plot(projected(1,i:nzcols:end), projected(2,i:nzcols:end));
subplot(3,3,4)
plot(projected(1:2,i:nzcols:end)');
subplot(3,3,5)
bar(selected(:,t));
ylim([0, 0.07])

idx2d = nzidx(i)-1;
x = mod(idx2d, length(Ws));
y = mod(idx2d/length(Ws), length(Hs));
z = floor(idx2d/length(Hs)/length(Ws));
idx3d = x+W*(y+H*z);
data=zeros(1,ts);
for i=1:ts+period-1
    filename = sprintf(strcat(filepath, filepattern), id_start+id_step*i )
    fp = fopen(filename, 'rb');
    fseek(fp, idx3d*4, 'bof');
    % a = fread(fp, [W*H*D,1], 'float32');
    % data(i) = a(idx3d+1);
    data(i) = fread(fp, 1, 'float32');
    fclose(fp);
end
subplot(3,3,6)
plot(data)
subplot(3,3,7)
plot(data(t:t+period-1));