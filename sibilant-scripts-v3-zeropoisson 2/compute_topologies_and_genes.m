%compute integrals and calculate topologies

%inputs:
%tfdata --> log-transformed data. rows = genes. columns = replicates (several per cell type). 
%idxs --> which cell types are we interested in out of all of the cell types in tfdata?
%labs --> cell type names
%iunique --> what is the cell type identity of each column of tfdata?
%iunique has same length as the # of columns in tfdata, and the max value
%of iunique is the number of cell types.

disp('Sibilant started.');

Params = get_parameters();
if Params.dist == 'P'
    if Params.sigmax > 1
        error('Params.sigmax for Zero Inflated Poisson distribution must be in the range [0,1]')
    end
end

[loggenemeans, loggenestds] = calc_log_mean_std(tfdata, iunique,Params.dist); 
%the function calc_log_mean_std has been modified to be able to calculate
%sigma and mu in ZIP distribution. However, these values will be used only
%if kde is on or filter gene by expression is on. 

disp('Beginning integration.');
[Integrals, best_mu, best_sigma,combinations] = calculate_integrals(tfdata, idxs, iunique, Params);
%combiantions contains all possible combinations of 3

ncells = length(idxs); ngenes = size(tfdata, 1);

%ncandid = # of likely topologies
%tlikely = most likely root of the triplet
%plikely = probability of T being the root

% Initialize holders
ncandid = zeros(size(combinations,1),1);
tlikely = zeros(size(ncandid));
plikely = zeros(size(ncandid));
pT_g_sum = zeros(size(ncandid,1), 5);

makeplot = 1; %make plots?
%find most likely topologies for each triplet

% Add some progress reporting
disp('Starting trajectories.');    
for j=1:size(combinations,1)
    icomb = j;
    iii = combinations(icomb,:);
    trip_data = loggenemeans(:,idxs(iii));
    [tlikely_ind, plikely_ind, ncandidates, pT_g_sum_j] = process_one_triplet(icomb, Integrals, Params, trip_data, makeplot, plot_dir, labs(iii));
    tlikely(j) = tlikely_ind;
    plikely(j) = plikely_ind;
    ncandid(j) = ncandidates;
    pT_g_sum(j,:) = pT_g_sum_j;
end

disp('Starting gene selection.');
%find union of marker and transition genes over all triplets
good_genes = zeros(ngenes,1);
transition_probs = zeros(ngenes, size(combinations,1));
marker_probs = zeros(ngenes, 3, size(combinations,1));
for i=1:size(combinations,1)
    if and(plikely(i)>Params.plikely_thresh,ncandid(i)==1)
        icomb = i;
        iii = combinations(icomb,:);
        trip_data = loggenemeans(:,idxs(iii));
        [prob_asym_topol, prob_marker_topol] = calculate_gene_probabilities(icomb,Integrals, Params,trip_data,tlikely(icomb));
        transition_probs(:,i) = prob_asym_topol;
        marker_probs(:,:,i) = prob_marker_topol;
        genestatus = find_classes(prob_asym_topol, prob_marker_topol, tlikely(icomb), trip_data, Params);
        good_genes = (good_genes + abs(genestatus))>0;
    end
end



