function Params = get_parameters()

%numeric integration parameters
Params.use_kde = 0; % do NOT turn this on until function is fixed to deal with n=[sbins,mbins]. 1 -> use kde-estimated empirical prior;  0 -> use uniform prior
Params.use_parallel_integrals = 1; %compute integrals in parallel?
Params.use_parallel_topology = 0;
Params.use_parallel_genes = 0;
Params.use_parallel = any([Params.use_parallel_integrals, Params.use_parallel_topology, Params.use_parallel_genes]);
Params.parallel_pool_size = 16; %number of processors for parallelization

Params.dist = 'N'; %distribution for gene expression. Default is 'N' for normal. Can also use 'P' for ZIP.
Params.mumin = 0; % min and max for mean prior in normal distribution, or for lambda prior in ZIP. 
Params.mumax = 7;
Params.sigmin = 0; % min and max for std prior in normal distribution, or for sigma prior in ZIP.
Params.sigmax = 3;
Params.mbins = 2^8; %number of bins for numeric integration over distribution parameters
Params.sbins = 2^8;

Params.cutoff = log2(10); % expression threshold for numeric integration

Params.verbose = 2; % 0 = default (no) reporting
                    % 1 = extra results
                    % 2 = debugging
Params.save_integrals = 1; % save integrals for later dot-plotting
Params.save_gene_probs = 1; % save gene probabilities for later
Params.save_best_dists = 1; % save info about 'best' mu and sigma per gene/cluster

Params.psla = [1/2 1/6 1/6 1/6]; %vector of prior probabilities for non-transition genes
Params.topology_prior = [1/4 1/4 1/4 1/4]; %topology prior [p(A) p(B) p(C) p(0)]
Params.plikely_thresh = 0.6; % threshold for probability of topology given data
Params.filter_min = 0; %ignore genes with values under threshold?
Params.log_asym_min = log2(2); %threshold for transition gene expression
Params.log_sym_min = log2(2); %threshold for marker gene expression
Params.oddsrange = logspace(-4,2,50); %range to vary prior odds p(beta_i = 1)/p(beta_i = 0)

%parameters for gene probabilities
Params.N_thresh_sym = 1000; % max number of genes in gene classes
Params.N_thresh_asym = 1000; 
Params.p_thresh_asym = 0.8; % probability thresholds for defining gene classes
Params.p_thresh_sym = 0.8; 
Params.odds0 = 0.05; % prior transition 

end