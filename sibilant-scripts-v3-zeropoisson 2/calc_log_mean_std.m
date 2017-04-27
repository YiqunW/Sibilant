function [loggenemeans, loggenestds] = calc_log_mean_std(data, iunique, g_dist)

loggenemeans = zeros(size(data,1), max(iunique));
loggenestds = zeros(size(data,1), max(iunique));

for i=1:max(iunique)
    datai= data(:,ismember(iunique,i));
    splmean = mean(datai,2);
    splstd = std(datai,[],2);
    if g_dist == 'N'
        loggenemeans(:,i) = splmean;
        loggenestds(:,i) = splstd;
    elseif g_dist == 'P'
        %maximum likelihood estimator of parameters:
        maxlbd = max(datai,[],2);
        d = (maxlbd-splmean)/500;
        lbd = repmat(d(:),1,500).*repmat(500:-1:1,size(data,1),1);
        lbd = lbd + repmat(splmean,1,500);
        F0= sum(datai==0,2)/size(datai,2);
        F0 = F0(:);
        a = (1-F0)./splmean(:);
        a = repmat(a(:),1,size(lbd,2));
        lhs = 1-a.*lbd;
        rhs = exp(-lbd);
        eq = lhs-rhs;
        [~,lbd0_ind] = min(abs(eq),[],2);
        lbd0 = diag(lbd(:,lbd0_ind));
        lbd0(F0==1) = 0;
        lbd0(a(:,1)~=a(:,1)) = 0;
        loggenemeans(:,i) = lbd0(:);
        loggenestds(:,i) = 1-splmean./lbd0(:);
        loggenestds(lbd0==0,i) = 1; %meaningless.
        
%         %method of moments estimator of parameters:
%         v1 = splstd.^2+splmean.^2-splmean;
%         v2 = splstd.^2-splmean;
%         disp([v1,v2])
%         loggenemeans(:,i) = v1./splmean;  
%         loggenestds(:,i) = v2./v1;
%         ind1 = find(splmean == 0);
%         loggenemeans(ind1,i) = 0;
%         loggenestds(ind1,i) = 1; %meaningless.
%         ind2 = find(v1<=0);
%         loggenemeans(ind2,i) = 0;
%         loggenestds(ind2,:) = 0;
%         loggenestds(v2<=0,i) = 0;
    end
end

end
