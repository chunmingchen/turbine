idx = 454348
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
x = mod(idx2d, length(Ws))+1;
y = mod(floor(idx2d/length(Ws)), length(Hs))+1;
z = floor(idx2d/length(Hs)/length(Ws))+1;
x3d = Ws(x)-1;
y3d = Hs(y)-1;
z3d = Ds(z)-1;
idx3d = x3d+W*(y3d+H*z3d);
data=zeros(1,ts);
for ii=1:ts+period-1
    filename = sprintf(strcat(filepath, filepattern), id_start+id_step*ii )
    fp = fopen(filename, 'rb');
    fseek(fp, idx3d*4, 'bof');
    % a = fread(fp, [W*H*D,1], 'float32');
    % data(i) = a(idx3d+1);
    data(ii) = fread(fp, 1, 'float32');
    fclose(fp);
end
subplot(3,3,6)
plot(data)
subplot(3,3,7)
plot(data(t:t+period-1));