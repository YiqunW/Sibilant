function [Integrals, best_mu, best_sigma,combinations] = calculate_integrals(log2tfdata, idxs, iunique, Params)

nbins = [Params.sbins, Params.mbins];
g_dist = Params.dist;
mumin = Params.mumin; mumax = Params.mumax;
sigmin = Params.sigmin; sigmax = Params.sigmax;
cutoff = Params.cutoff;
use_parallel_integrals = Params.use_parallel_integrals;
[loggenemeans, loggenestds] = calc_log_mean_std(log2tfdata, iunique, g_dist);

%estimate mean and standard deviation prior
if Params.use_kde == 1 %Estimate empirical distribution of means and standard deviations
    [~,density,M,S]=kde2d([loggenemeans(:) loggenestds(:)],nbins,[mumin sigmin],[mumax sigmax]); %kernel density estimation.
    % Reference: Botev. Z.I., Grotowski J.F and Kroese D. P. (2010). Kernel density estimation via diffusion. Annals of Statistics.
    % Volume 38, Number 5, Pages 2916--2957
    density(density(:)<0) = 0; %get rid of numerical errors (sum of negative elements was ~10^-8 in testing)
    dm=M(1,2)-M(1,1);
    ds=S(2,1)-S(1,1);
    P = density/(dm*ds*sum(density(:)));
else %uniform prior over means and standard deviations
    dm = (mumax-mumin)/nbins(2);
    ds = (sigmax-sigmin)/nbins(1); % will be zero is sigmax==sigmin
    if sigmax==sigmin
        ds=1;
    end
    [M,S]=meshgrid(mumin:dm:mumax,sigmin:ds:sigmax);
    density = ones(size(M));
    P = density/(dm*ds*sum(density(:)));
end


% Compute integrals

ngenes = size(log2tfdata,1);
ncells = length(idxs);
combinations = combnk(1:ncells,3);
lM = size(M,1);
lS = size(S,2);
Mu = M(1,2:end);
firstcut = find(Mu>cutoff,1,'first');

IA = zeros(ngenes, ncells);%integrals over one cell type (to use in integrals for Equation 8)
IAB = zeros(ncells,ncells,ngenes);%integrals over 2 cell types (to use in integrals for Equation 8)
IABCsame = zeros(ngenes, size(combinations,1));%integrals over 3 cell types (Equation 9) muA=muB=muC
IA_BC_Amax = zeros(ngenes, size(combinations,1));%integrals with muA < muB=muC (Equation 8)
IA_BC_Bmax = zeros(ngenes, size(combinations,1));%integrals with muB < muA=muC (Equation 8)
IA_BC_Cmax = zeros(ngenes, size(combinations,1));%integrals with muC < muB=muA (Equation 8)
IABC_Amin = zeros(ngenes, size(combinations,1));%integrals with muA < muB and muA < muC (Equation 5)
IABC_Bmin = zeros(ngenes, size(combinations,1));%integrals with muB minimum (Equation 5)
IABC_Cmin = zeros(ngenes, size(combinations,1));%integrals with muC minimum (Equation 5)

best_mu = zeros(ngenes, ncells);
best_sigma = zeros(ngenes, ncells);

dm=M(1,2)-M(1,1);
ds=S(2,1)-S(1,1);

if Params.use_parallel
    %start parallel cluster
    pc = parcluster('local');
    %pc.JobStorageLocation = strcat('/scratch/furchtg/',getenv('SLURM_JOB_ID'));
    parpool(pc,Params.parallel_pool_size)
end

%% if assume normal distribution
if g_dist=='N'
    logS = log(2*pi*S.^2);
elseif g_dist=='P'
    logS1 = log(1-S);
    ExpM = exp(-M); %exp(-mu). sbins by mbins matrix.
    logSM = log(S+(1-S).*ExpM); %log[sigma+(1-sigma)exp(-mu)]. sbins by mbins matrix.
    logS = {logS1,logSM};
else
    disp('Invalid distribution parameter. Acceptable input: ''N'' or ''P''')
end

if use_parallel_integrals
    parfor i=1:ngenes
        if (mod(i,25) == 0)
            tic;
        end
        [tempIA, tempIAB, tempIABCsame, tempIABC_Amin, tempIABC_Bmin, tempIABC_Cmin, tempIA_BC_Amax, tempIA_BC_Bmax, tempIA_BC_Cmax, best_mu_gene, best_sigma_gene] = calculate_integrals_one_gene(g_dist, i, ncells, combinations, log2tfdata, iunique, idxs, P, M, S, logS, lM, lS, firstcut, dm, ds);
        IA(i,:) = tempIA;
        IAB(:,:,i) = tempIAB;
        IABCsame(i,:) = tempIABCsame;
        IABC_Amin(i,:) = tempIABC_Amin;
        IABC_Bmin(i,:) = tempIABC_Bmin;
        IABC_Cmin(i,:) = tempIABC_Cmin;
        IA_BC_Amax(i,:) = tempIA_BC_Amax;
        IA_BC_Bmax(i,:) = tempIA_BC_Bmax;
        IA_BC_Cmax(i,:) = tempIA_BC_Cmax;
        best_mu(i,:) = best_mu_gene;
        best_sigma(i,:) = best_sigma_gene;
        if (mod(i,25) == 0)
            disp(strcat(num2str(i), 'th gene integrated and took ', num2str(toc), ' seconds.'));
        end
    end
    clear pc
else
    for i=1:ngenes
        [tempIA, tempIAB, tempIABCsame, tempIABC_Amin, tempIABC_Bmin, tempIABC_Cmin, tempIA_BC_Amax, tempIA_BC_Bmax, tempIA_BC_Cmax, best_mu_gene, best_sigma_gene] = calculate_integrals_one_gene(g_dist, i, ncells, combinations, log2tfdata, iunique, idxs, P, M, S, logS, lM, lS, firstcut, dm, ds);
        IA(i,:) = tempIA;
        IAB(:,:,i) = tempIAB;
        IABCsame(i,:) = tempIABCsame;
        IABC_Amin(i,:) = tempIABC_Amin;
        IABC_Bmin(i,:) = tempIABC_Bmin;
        IABC_Cmin(i,:) = tempIABC_Cmin;
        IA_BC_Amax(i,:) = tempIA_BC_Amax;
        IA_BC_Bmax(i,:) = tempIA_BC_Bmax;
        IA_BC_Cmax(i,:) = tempIA_BC_Cmax;
        best_mu(i,:) = best_mu_gene;
        best_sigma(i,:) = best_sigma_gene;
        if (mod(i,100) == 0)
            disp(strcat(num2str(i), ' genes integrated.'));
        end
    end
end

Integrals.IA = IA;
Integrals.IAB = IAB;
Integrals.IABCsame = IABCsame;
Integrals.IABC_Amin = IABC_Amin;
Integrals.IABC_Bmin = IABC_Bmin;
Integrals.IABC_Cmin = IABC_Cmin;
Integrals.IA_BC_Amax = IA_BC_Amax;
Integrals.IA_BC_Bmax = IA_BC_Bmax;
Integrals.IA_BC_Cmax = IA_BC_Cmax;

end
function F1 = logDa(g_dist,M,S,xk1,logS,lM,lS)
%%This function returns a sbins by mbins matrix of logDa(gi|mu,sigma) for gene(i) in cluster j.
%input arguments:
% g_dist: Assumption for gene expression.'N' for normal, or 'P' for ZIP (zero inflated poisson).
% M,S: meshgrid components (each is a sbin by mbin matrix) created with the prior mu and sigma values.
% xk1: gene(i) expression levels in all cells in cluster j
% logS: Terms independent of gi in logDa: Equals log[2pi(sigma)^2] for
% 'N', and {log(1-sigma), log[sigma+(1-sigma)exp(-mu)]} for 'P'.
% lM, lS: number of prior mu and sigma values. Should be roughly mbins, sbins
N1 = length(xk1); %number of cells in cluster j
if g_dist == 'N'
    F1 = -N1*logS/2; %-(N/2)log[2pi(sigma)^2], sbins by mbins matrix, same value for each row.
    [MM,XK] = meshgrid(M(1,:),xk1); %MM, XK are Ncells by mbins matrices. MM contains prior mu values column-wise. XK contains gene(i) levels row-wise.
    F1a = repmat(sum(-(MM-XK).^2,1),[lM,1]); %sum[-(gi-mu)^2]. sbins by mbins matrix. Identical rows.
    F1 = F1 + F1a./(2*S.^2); %logDa = -(N/2)log[2pi(sigma)^2]+sum[-(gi-mu)^2]/(2sigma^2). sbins by mbins.
elseif g_dist == 'P'
    logS1 = logS{1}; %log(1-sigma). sbins by mbins matrix. Identical rows.
    logSM = logS{2}; %log[sigma+(1-sigma)exp(-mu)]. sbins by mbins matrix.
    xk1p = xk1(xk1>0); %gene(i) levels in all cells expression gene(i) in cluster j. Vector.
    Np = length(xk1p); %number of cells with positive gene(i) expression.
    N0 = N1-Np; %number of cells with zero gene(i) level.
    sum_g = sum(xk1p); %single value
    sum_gamma = sum(gammaln(xk1p+1)); %for gi>0, sum{log[gamma(gi+1)]}. Single value.
    sum_Mg = log(M(1,:))*sum_g; %vector of length mbins. log(mu)*Sum(gi) for gi>0.
    F1 = Np*logS1+N0*logSM-Np*M; %Np*log(1-sigma)+N0*log[sigma+(1-sigma)exp(-mu)]-Np*mu. sbins by mbins matrix.
    F1 = F1 + repmat(sum_Mg, [lS,1]); %F1 as above + log[mu*Sum(gi) for gi>0]. sbins by mbins matrix.
    F1 = F1 - sum_gamma; %logDa = Np*log(1-sigma)+N0*log[sigma+(1-sigma)exp(-mu)]-Np*mu+log[mu*Sum(gi)]-sum{log[gamma(gi+1)]}, for gi>0.
end
end
function [tempIA, tempIAB, tempIABCsame, tempIABC_Amin, tempIABC_Bmin, tempIABC_Cmin, tempIA_BC_Amax, tempIA_BC_Bmax, tempIA_BC_Cmax, best_mu_gene, best_sigma_gene] = calculate_integrals_one_gene(g_dist, i, ncells, combinations, log2tfdata, iunique, idxs, P, M, S, logS, lM, lS, firstcut, dm, ds)
tempIA = zeros(1,ncells);
tempIAB = zeros(ncells,ncells,1);
tempIABCsame = zeros(1, size(combinations,1));
tempIA_BC_Amax = zeros(1, size(combinations,1));
tempIA_BC_Bmax = zeros(1, size(combinations,1));
tempIA_BC_Cmax = zeros(1, size(combinations,1));
tempIABC_Amin = zeros(1, size(combinations,1));
tempIABC_Bmin = zeros(1, size(combinations,1));
tempIABC_Cmin = zeros(1, size(combinations,1));
Isigma = cell(ncells,1); %pre-allocate vector for integration over sigma for one cell type
Isigma2 = cell(ncells,ncells); %pre-allocate matrix for integration over sigma for two cell types
diags = cell(ncells,ncells); %pre-allocate cell for integration over mu for two cell types
data = log2tfdata(i,:); %vector of gene(i) expression levels in all cells
best_mu_gene = zeros(1,ncells); %pre-allocate vector for the parameters that give the highest likely hood for each gene in each cell type
best_sigma_gene = zeros(1,ncells);
for j=1:ncells %for each cell cluster j
    %N1 = sum(iunique==idxs(j));%integrals over one cell type (to use in integrals for Equation 9)
    xk1 = data(ismember(iunique,idxs(j))); %gene(i) expression levels in all cells in cluster j
    F1 = logDa(g_dist,M,S,xk1,logS,lM,lS); %returns a sbin by mbin matrix of logDa(gi|mu,sigma) for gene(i) in cluster j. See equation (11) (or the equivalent for ZIP).
    F1 = exp(F1);
    F1 = P.*F1;  % multiply by prior, but flat.
    F1(F1~=F1) = 0; % correct for NaN?
    
    % Save most likely distribution for each gene, per cluster.
    [ind_r, ind_c]=find(F1==max(F1(:)));
    best_mu_gene(j) = M(1,ind_c(1));
    best_sigma_gene(j) = S(ind_r(1),1);
    
    tempIA(1,j) = dm*ds*trapz(trapz(F1,1),2); %integrate Da over all prior mu and sigma for gene(i) in cluster j. Can be used in Eq.7 for T=0.
    Isigma{j} = ds*trapz(F1,1);%store integrals over sigma for use later. Isigma{j} is a vector of length mbins.
    for k=(j+1):ncells %integrals over 2 cell types (to use in integrals for Equation 9)
        %N2 = sum(iunique==idxs(k)); N12 = N1+N2;
        xk2 = data(ismember(iunique,idxs(k)));
        xk12 = [xk1 xk2];
        F2 = logDa(g_dist,M,S,xk12,logS,lM,lS);
        F2 = exp(F2);
        F2(F2~=F2) = 0;
        F2 = P.*F2;
        tempIAB(j,k,1) = dm*ds*trapz(trapz(F2,1),2); %integrate Da over the hyperprior space for gene(i) in cluster j+k.
        tempIAB(k,j,1) = tempIAB(j,k,1); %tempIAB is a symmetric 3by3 matrix.
        Isigma2{j,k} = ds*trapz(F2,1); %Isigma2 is a symmetric 3by3 cell.
        Isigma2{k,j} = Isigma2{j,k};
    end
end
for j=1:ncells
    for k=(j+1):ncells
        AB = repmat(reshape(Isigma{j},[lM 1]),[1 lM]).*repmat(reshape(Isigma{k},[1 lM]),[lM 1]);
        AB_diag = diag(cumtrapz(cumtrapz(AB(lM:-1:1,lM:-1:1),1),2));
        AB_diag = AB_diag(lM:-1:1); %AB_diag is a vector of length mbins.
        %Each element is the integral of [Integral(Dj over sigma)*Integral(Dk over sigma)]
        %over mu values greater than a certain mu (in the prior space).
        diagjk = AB_diag(2:end)'; %delete the first element of integration of Dj and Dk over the complete hyperprior space.
        diagjk(1:(firstcut-1)) = diagjk(firstcut); %replace integrations over mu>mu_replace for mu_replace<mu_cutoff with the integrations over mu>mu_cutoff
        diags{j,k} = diagjk; %diags is a 3by3 symmetric cell with each element being a (mbins-1)by1 vector
        diags{k,j} = diags{j,k};
    end
end

for j=1:size(combinations,1)
    xk123 = data(ismember(iunique,idxs(combinations(j,:)))); %computation of integrals over 3 cell types assuming they are from a single distribution
    %N123 = length(xk123);
    F3 = logDa(g_dist, M,S,xk123,logS,lM,lS);
    F3 = exp(F3);
    F3 = P.*F3;
    F3(F3~=F3) = 0;
    tempIABCsame(1,j) = dm*ds*trapz(trapz(F3,1),2); % Equation (10). Integration of Da over the entire hyperprior space
    
    tempIABC_Amin(1,j) = dm^3*trapz(Isigma{combinations(j,1)}(1:lM-1).*diags{combinations(j,2),combinations(j,3)}); % Equation (6)
    tempIABC_Bmin(1,j) = dm^3*trapz(Isigma{combinations(j,2)}(1:lM-1).*diags{combinations(j,1),combinations(j,3)}); % if cutoff=0, this is exactly Equation(6)
    tempIABC_Cmin(1,j) = dm^3*trapz(Isigma{combinations(j,3)}(1:lM-1).*diags{combinations(j,1),combinations(j,2)});
    
    
    A_BC = repmat(reshape(Isigma{combinations(j,1)},[lM 1]),[1 lM]).*repmat(reshape(Isigma2{combinations(j,2),combinations(j,3)},[1 lM]),[lM 1]);
    B_AC = repmat(reshape(Isigma{combinations(j,2)},[lM 1]),[1 lM]).*repmat(reshape(Isigma2{combinations(j,1),combinations(j,3)},[1 lM]),[lM 1]);
    C_AB = repmat(reshape(Isigma{combinations(j,3)},[lM 1]),[1 lM]).*repmat(reshape(Isigma2{combinations(j,1),combinations(j,2)},[1 lM]),[lM 1]);
    A_BC = tril(A_BC,1);
    B_AC = tril(B_AC,1);
    C_AB = tril(C_AB,1);
    tempIA_BC_Amax(1,j) = dm^2*trapz(trapz(A_BC,1),2); % Equation (9) but why tril(X,1) instead of tril(X,0)??
    tempIA_BC_Bmax(1,j) = dm^2*trapz(trapz(B_AC,1),2);
    tempIA_BC_Cmax(1,j) = dm^2*trapz(trapz(C_AB,1),2);
end
end