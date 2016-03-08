RUN_GRUBBS_TEST = 0
RUN_Mahalanobis = 0
RUN_KLDivergence = 0
RUN_EMD = 1
RUN_MEAN = 0

USE_GMM = 1

CASE=1 
if CASE == 1
%     filepath = '/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/';      id_start = 49401
%     filepath = '/home/chenchu/volumes/TEMP/13.80_141219/';     id_start = 6201
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00/';    id_start = 20601
    
    filepattern = 's35_noinj.r2b%d.p3d.g%d_Pressure.raw'
    output_filepattern_grubbs = 'distribution_diff/distdiff_b%d_%d.raw'
    output_filepattern_mahalanobis = 'distribution_diff/mahalanobis_b%d_%d.raw'
    output_filepattern_KLD = 'distribution_diff/kldiv_b%d_%d.raw'
    output_filepattern_EMD = 'distribution_diff/emd_b%d_%d.raw'
    output_filepattern_mean = 'distribution_diff/mean_b%d_%d.raw'
    id_step = 25
    nfiles = 576
    W=151;  % Note!! Original range
    H=71
    D=56;
    res = 5;
end

NBINS = 64;
nblocks = 36;
maxval = 4;
edges = linspace(0, maxval, NBINS+1);
w1 = floor(W/res);
h1 = floor(H/res);
d1 = floor(D/res);
ncells = w1*h1*d1
Z_TH = 3.3296;
% center of each bin:
bin_val = (edges(1:end-1)+edges(2:end)) * 0.5; bin_val = bin_val';
% preallocate necessary buffers
histmat = zeros(NBINS, ncells, nblocks);
anomaly = zeros(ncells,nblocks);
anomaly_count = zeros(ncells,nblocks);
gmmmat = cell(ncells, nblocks);
for t=1:nfiles
    for b=1:nblocks
        filename = sprintf(strcat(filepath, filepattern), b, id_start+id_step*(t-1) )
        fp = fopen(filename, 'rb');
        data = fread(fp, W*H*D, 'float32');    
        fclose(fp);
        data = reshape(data, [W,H,D]);
        %% compute histogram
        count = 1;
        for z=1:res:d1*res
            for y=1:res:h1*res
                for x=1:res:w1*res
                    tmp = data(x:x+res-1, ...
                               y:y+res-1, ...
                               z:z+res-1 );
                    tmp = reshape(tmp, [res*res*res,1]);
                    if USE_GMM
                        gmm = isp_gmmtrain(tmp', 'nMixtures', 4, 'silent', true);
                        gmmmat{count,b} = gmm;
                        
%                         
%                         gmmmat{count,b} = gmdistribution.fit(tmp, 4,  'CovType','diagonal');
% %                         verfication:
%                         h = histc(tmp, edges);                     
%                         hold off
%                         plot( h )
%                         hold on
%                         plot( pdf(gmmmat{count,b}, edges') );
%                         cov = zeros(1,1,4);
%                         cov(1,1,:) = gmm.covariance;
%                         plot( pdf(gmdistribution(gmm.mean', cov, gmm.weight'), edges'))
%                         pause
%                         hold off
                        h = histc(tmp, edges);
                        histmat(:,count,b) = h(1:NBINS);
                    else
%                     hh = histcounts(tmp, edges);
%                     histmat(:,count,b) = hh';
                        h = histc(tmp, edges);
                        histmat(:,count,b) = h(1:NBINS);
                    end
                                        
                    count = count +1;
                end
            end
        end
    end
    
    %% anomaly analysis
    if RUN_GRUBBS_TEST 
        for i=1:ncells
            cell_hist = sum( histmat(:,i,:), 3 );
            npoints = sum(cell_hist);
            cell_mean = sum(bin_val.*cell_hist)/npoints;
            cell_std = sqrt(sum(bin_val.*bin_val.*cell_hist)/npoints - cell_mean*cell_mean); % sqrt(E(X^2)-E(X)^2)

            anomaly_bins = (bin_val > cell_mean + cell_std * Z_TH + maxval/NBINS*.5) |  ...
                           (bin_val < cell_mean - cell_std * Z_TH - maxval/NBINS*.5);
            for b=1:nblocks
                anomaly_count(i, b) = sum(anomaly_bins .* histmat(:,i,b) );
    %             if (anomaly_count(i,b) > 0)
    %                 anomaly_count(i,b)
    % %                 pause
    %             end 
            end        
        end
        for b=1:nblocks
            filename = strcat(filepath, sprintf(output_filepattern_grubbs, b, t) )
            fp = fopen(filename, 'wb')
            fwrite(fp, anomaly_count(:,b), 'float32');
            fclose(fp);
        end
    end
    % regard each distribution as a high-dim vector
    if RUN_Mahalanobis 
        for i=1:ncells
            mat = squeeze(histmat(:,i,:));  % [bins blocks]
            c = mean(mat, 2);
            matc = mat - repmat(c, 1, nblocks);
            
            % nonzero rows
            rows = find(all(matc,2));
            % all 0?
            if (length(rows) == 0)
                anomaly(i,:) = zeros(1, nblocks);
                continue;
            end
            matc1 = matc(rows, :);
            
            cov = matc1*matc1'/nblocks;
            
            % svd
            if 0
                [U,S,V] = svd(cov);
                s=diag(S);            
                s(find(s==0)) = 1;   
                for b=1:nblocks
                    offset = matc1(:,b);
                    anomaly(i,b) = sqrt( offset'*V*(U'*offset./s)  );
                end
            end
            % eigen analysis.  faster than svd
            if 1
                [PC, V] = eig(cov);
                V=diag(V);            
                V(find(V==0)) = 1;                      
                for b=1:nblocks
                    offset = matc1(:,b);
                    anomaly(i,b) = sqrt( sum((PC' * offset).^2 ./ V ) );
                end
            end
            % original definition. not stable
            if 0
                invcov = inv(cov);
                a= zeros(1,b);
                for b=1:nblocks
                    offset = matc1(:,b);
                    a(b) = sqrt(offset'*invcov*offset);
                end
                a
            end
        end
        
        for b=1:nblocks
            filename = strcat(filepath, sprintf(output_filepattern_mahalanobis, b, t) )
            fp = fopen(filename, 'wb')
            fwrite(fp, anomaly(:,b), 'float32');
            fclose(fp);
        end
    end
    
    % regard each distribution as a high-dim vector, and compare with 
    if RUN_KLDivergence
        for i=1:ncells
            mat = squeeze(histmat(:,i,:));  % [bins blocks]
            mat = mat/(res*res*res);        % get probability
            
            % mean distribution            
            q = mean(mat, 2); 
            q(q==0) = 1; % esacpe 0
            qmat = repmat(q, 1, nblocks);
            
            mat1 = mat;  mat1(mat1==0)=1;   % escape 0            
            d = sum(mat .* log(mat1 ./ qmat), 1);
            
            % imagesc([mat q; d 1])
            anomaly(i,:) = d;
            
        end
        
        for b=1:nblocks
            filename = strcat(filepath, sprintf(output_filepattern_KLD, b, t) )
            fp = fopen(filename, 'wb');
            fwrite(fp, anomaly(:,b), 'float32');
            fclose(fp);
        end
    end
    
    if RUN_EMD
        for i=1:ncells
            if ~USE_GMM
                mat = squeeze(histmat(:,i,:));  % [bins blocks]
                mat = mat/(res*res*res);        % get probability
                cummat = cumsum(mat, 1); % CDF

                % mean distribution            
                c = mean(mat, 2); 
                cumc = cumsum(c);
                cumcmat = repmat(cumc, 1, nblocks);

                anomaly(i, :) = sum( abs(cumcmat-cummat) );  % L1
%                 anomaly(i, :) = sqrt(sum( (cumcmat-cummat).^2 )); % L2
%             imagesc([mat c; anomaly(i,:) 0]);
%             pause
            else                
                means = zeros(1,36*4);
                vars = zeros(1,36*4);
                weights = zeros(1,36*4);
                for b=1:36
                    gmm = gmmmat{i,b}
                    means((b-1)*4+1 : (b-1)*4+4) = gmm.mean;
                    vars((b-1)*4+1 : (b-1)*4+4) = gmm.covariance;                    
                    weights((b-1)*4+1 : (b-1)*4+4) = gmm.weight;
                end
                gmm.mean = means;
                gmm.covariance = vars;
                gmm.weight = weights/36;
                gmm.nMixtures = 4*36;
                for b=1:36
                    d = isp_gmmdistance(gmm, gmmmat{i,b})
                    anomaly(i,b)=d;
                end
                
                % verify
                mat = squeeze(histmat(:,i,:));  % [bins blocks]
                mat = mat/(res*res*res);        % get probability         
                c = mean(mat, 2); 
                cummat = cumsum(mat, 1); % CDF                
                imagesc([mat c; anomaly(i,:) 0]);
                pause
            end
            
        end
        
        for b=1:nblocks
            filename = strcat(filepath, sprintf(output_filepattern_EMD, b, t) )
            fp = fopen(filename, 'wb')
            fwrite(fp, anomaly(:,b), 'float32');
            fclose(fp);
        end
    end
    
    if RUN_MEAN
        pts = res*res*res;
        for i=1:ncells
            mat = squeeze(histmat(:,i,:));  % [bins blocks]
            dist_average = sum(mat .* repmat(bin_val, 1, nblocks), 1) / pts;
            
            s = std(dist_average);
            m = mean(dist_average);
            if s==0
                anomaly(i,:) = zeros(1, nblocks);
                continue;
            end
            anomaly(i,:) = (dist_average' - m)/s;            
        end
        
        for b=1:nblocks
            filename = strcat(filepath, sprintf(output_filepattern_mean, b, t) )
            fp = fopen(filename, 'wb')
            fwrite(fp, anomaly(:,b), 'float32');
            fclose(fp);
        end
    end
    
end

