filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
filepattern = '14.2-13.8-turbine_%d.vti_Pressure_slice44.raw'
nfiles = 576
period = 144
W=509
H=509

X = 88
Y = 400
Z = 44

data = zeros(nfiles, 1);
for j=1:nfiles
    filename = sprintf(strcat(filepath, filepattern), j)

    fp = fopen(filename, 'rb');
    raw = fread(fp, [509,509], 'float32');
    data(j) = raw(Y,X);
    fclose(fp);
end

plot(data)