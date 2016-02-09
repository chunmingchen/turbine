function output = mycwt(data, widths)

output = zeros(length(widths), length(data));
for count = 1:length(widths)
    width = widths(count);
    wavelet_data = mymorlet(min(10*width, length(data)), width, 1);
    output(count, :) = conv(data, wavelet_data, 'same');
end

end