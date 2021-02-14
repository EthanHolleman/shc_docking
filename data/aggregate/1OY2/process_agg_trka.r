# Process results of a multi iteration run into one large score file
# that has paths to pdb files, distances to PTB domain and
# average ligand position. This aggregated dataframe is saved
# as an RDS file and can be downloaded to local machine for easier plotting
# (interactively) using Rstudio.

source('/home/ethollem/software/RDBC/handler/R/energy_bowl.r')
agg_results_file = '/home/ethollem/jobs/random_docking_trka/aggregate/1OY2/1OY2_agg.tsv'
results_dir = '/home/ethollem/jobs/random_docking_trka/results_1OY2'
output_file = '/home/ethollem/jobs/random_docking_trka/aggregate/1OY2/1OY2_ag.rds'
number_threads = 4

cmd <- 'python3 ~/software/RDBC/rh.py -o %s -mia %s'
cmd <- sprintf(cmd, results_dir, agg_results_file)
print(cmd)
system(cmd)
print('Reading score file')
score_file.df <- read_tab_score_file(agg_results_file)
print('Removing uncontacted ligands')
score_file.df <- remove_uncontacted_ligands(score_file.df)
print('Adding pdb paths')
score_file.df <- add_pdb_path_for_multi_iter_run(
    score_file.df, results_dir
)
saveRDS(score_file.df, output_file)
print('Calculting average chain coords')
score_file.df <- add_average_chain_coords_multi_iter(
    score_file.df, threads=number_threads)
saveRDS(score_file.df, output_file)
score_file.df <- distance_to_binding_pocket(score_file.df)
saveRDS(score_file.df, output_file)

