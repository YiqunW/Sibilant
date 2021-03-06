% Must pass sample_dir in the calling script
% Hard-coded directory parameters
base_dir = pwd;
sample_dir = '/simuN';
script_dir = '/';

plot_dir = strcat(base_dir, sample_dir, '/plots/');
disp(strcat('Loading data: ', sample_dir));

% Load data
cd(strcat(base_dir, sample_dir));
mkdir('plots');
tfdata = dlmread('tfdata.txt');
idxs = importdata('idxs.txt');
labs = importdata('labs.txt');
if ~iscell(labs)
    labs = cellfun(@num2str, num2cell(labs), 'UniformOutput', 0);
end
iunique = importdata('iunique.txt');

% Run sibilant
disp('Starting sibilant.');
cd(strcat(base_dir, script_dir));
run('compute_topologies_and_genes.m');

% Write out the data it produces
disp('Saving data.');
cd(strcat(base_dir, sample_dir));
topologies = table(combinations, pT_g_sum, ncandid, tlikely, plikely);
writetable(topologies,'topologies.txt','Delimiter','\t');
dlmwrite('goodgenes.txt',good_genes,'\n');

% Write out optional data if desired
disp('Saving optional data.');
if Params.save_integrals
    save('integrals.mat', 'Integrals');
end
if Params.save_gene_probs
    save('geneprobs_marker.mat', 'marker_probs');
    save('geneprobs_transition.mat', 'transition_probs');
end
if Params.save_best_dists
    save('bestdist_mu.mat', 'best_mu');
    save('bestdist_sigma.mat', 'best_sigma');
end
disp('All done.');