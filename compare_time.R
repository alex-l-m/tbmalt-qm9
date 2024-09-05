# Compare runtime of the three methods used: DFTB+, tbmalt individual, tbmalt
# batched, and make histograms
library(tidyverse)
library(cowplot)
library(ggdark)
library(glue)
theme_set(dark_mode(theme_cowplot(font_size = 24)))

tbmalt_result_cols <- cols(
    molecule = col_character(),
    status = col_character(),
    energy = col_double(),
    repulsive_energy = col_double(),
    scc_energy = col_double(),
    run_time = col_double()
)

tbmalt_results_individually <- read_csv('tbmalt_results_individually.csv', col_types = tbmalt_result_cols) |>
    rename(mol_id = molecule)

dftbplus_results <- read_csv('ase_dftb_qm9.csv', col_types = cols(
    mol_id = col_character(),
    energy = col_double(),
    run_time = col_double()
))

run_times <- tbmalt_results_individually |>
    select(mol_id, run_time) |>
    inner_join(dftbplus_results, by = 'mol_id', suffix = c('_tbmalt', '_dftb')) |>
    pivot_longer(cols = c(run_time_tbmalt, run_time_dftb), names_to = 'method', values_to = 'run_time') |>
    mutate(method = str_replace(method, 'run_time_', ''))


histograms <- run_times |>
    ggplot(aes(x = run_time)) +
    facet_wrap(~method, ncol = 1) +
    geom_histogram(binwidth = .01)

ggsave('histograms.png', histograms, height = unit(11.51, 'in'), width = unit(4.75, 'in'))

tbmalt_results_batch <- read_csv('tbmalt_results.csv', col_types = tbmalt_result_cols) |>
    rename(mol_id = molecule)
run_times_incl_batch <- bind_rows(
    tbmalt_individually = filter(tbmalt_results_individually, status == 'success'),
    dftbplus = dftbplus_results,
    tbmalt_batch = tbmalt_results_batch,
    .id = 'method') |>
    # Include only molecule id's present in all three
    filter(mol_id %in% dftbplus_results$mol_id &
           mol_id %in% filter(tbmalt_results_individually, status == 'success')$mol_id &
           mol_id %in% tbmalt_results_batch$mol_id)
avg_run_times_tbl <- run_times_incl_batch |>
    group_by(method) |>
    summarize(avg_run_time = mean(run_time),
              total_run_time = sum(run_time),
              n_mols = length(run_time),
              .groups = 'drop')
print(avg_run_times_tbl)
