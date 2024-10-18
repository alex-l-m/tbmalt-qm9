library(tidyverse)
library(robustbase)
library(cowplot)
library(ggdark)
library(tidymodels)
library(glue)
theme_set(dark_mode(theme_cowplot(font_size = 16)))

dftb_results_infile <- 'ase_dftb_qm9.csv'

dftb <- read_csv(dftb_results_infile, col_types = 
cols(
  mol_id = col_character(),
  energy = col_double())) |>
    rename(dftb_energy = energy)

training_output <- read_csv('trained.csv', col_types = cols(
    epoch = col_integer(),
    mol_id = col_character(),
    energy = col_double(),
    repulsive_energy = col_double(),
    scc_energy = col_double(),
    run_time = col_double(),
    split = col_character()
))

trained_results <- training_output |>
    filter(epoch == max(epoch, na.rm = TRUE)) |>
    rename(trained_energy = energy) |>
    select(mol_id, trained_energy, split)

untrained_results <- training_output |>
    filter(epoch == 0) |>
    rename(tbmalt_energy = energy) |>
    select(mol_id, tbmalt_energy)

TO_EV = 27.211386246
pre_regression_table <- read_csv('regression_table.csv', col_types = cols(
    mol_id = col_character(),
    qm9_id = col_integer(),
    n_H = col_double(),
    n_C = col_double(),
    n_N = col_double(),
    n_O = col_double(),
    qm9_u0 = col_double()
)) |>
    # Trained tbmalt
    inner_join(trained_results, by = 'mol_id') |>
    # Untrained tbmalt
    inner_join(untrained_results, by = 'mol_id') |>
    # Convert electronvolts
    mutate(qm9_u0 = qm9_u0 * TO_EV,
           trained_energy = trained_energy * TO_EV,
           tbmalt_energy = tbmalt_energy * TO_EV) |>
    # Total number of atoms
    mutate(n_atoms = n_H + n_C + n_N + n_O) |>
    # Join with DFTB results
    inner_join(dftb, by = 'mol_id') |>
    # Filter out NA rows
    filter(!is.na(tbmalt_energy), !is.na(qm9_u0), !is.na(dftb_energy), !is.na(trained_energy)) |>
    # Subtract out atomic contributions with regression
    mutate(tbmalt_energy_ref = predict(lmrob(tbmalt_energy ~ n_H + n_C + n_N + n_O)),
           qm9_u0_ref = predict(lmrob(qm9_u0 ~ n_H + n_C + n_N + n_O)),
           dftb_energy_ref = predict(lmrob(dftb_energy ~ n_H + n_C + n_N + n_O)),
           trained_energy_ref = predict(lmrob(trained_energy ~ n_H + n_C + n_N + n_O)),
           trained_energy_ref = tbmalt_energy_ref,
           tbmalt_energy_diff = tbmalt_energy - tbmalt_energy_ref,
           qm9_u0_diff = qm9_u0 - qm9_u0_ref,
           dftb_energy_diff = dftb_energy - dftb_energy_ref,
           trained_energy_diff = trained_energy - trained_energy_ref)

# Calculate the range of the energy per atom in QM9
lower_limit <- with(pre_regression_table, min(qm9_u0_diff / n_atoms))
upper_limit <- with(pre_regression_table, max(qm9_u0_diff / n_atoms))

regression_plot_tbl <- pre_regression_table |>
    filter(tbmalt_energy_diff / n_atoms > lower_limit - .2) |>
    filter(tbmalt_energy_diff / n_atoms < upper_limit + .2)

regression_plot <- regression_plot_tbl |>
    ggplot(aes(y = tbmalt_energy_diff / n_atoms, x = qm9_u0_diff / n_atoms)) +
    # y=x line
    geom_abline(intercept = 0, slope = 1, color = 'blue', linetype = 'dashed') +
    geom_point(size = 0, alpha = 0.5) +
    geom_smooth(method = lmrob, se = FALSE) + 
    coord_obs_pred() +
    ylab('tbmalt (not trained)') + xlab('QM9 DFT')

pass_filter <- nrow(regression_plot_tbl)
all <- nrow(pre_regression_table)
print(glue('Fraction not shown: {1 - pass_filter / all}'))

combined_regression_tbl <- pre_regression_table |>
    # Normalize by number of atoms
    mutate(`Untrained TBMaLT` = tbmalt_energy_diff / n_atoms,
           qm9_u0 = qm9_u0_diff / n_atoms,
           `Standard DFTB (DFTB+)` = dftb_energy_diff / n_atoms,
           `Trained TBMaLT` = trained_energy_diff / n_atoms) |>
    pivot_longer(cols = c(`Untrained TBMaLT`, `Standard DFTB (DFTB+)`, `Trained TBMaLT`),
                 names_to = 'method', values_to = 'energy_diff_peratom') |>
    # Put standard first in the plot, by converting to a factor
    mutate(method = factor(method, levels = c('Standard DFTB (DFTB+)',
                                              'Untrained TBMaLT',
                                              'Trained TBMaLT'))) |>
    mutate(outlier = energy_diff_peratom - qm9_u0 > 0.2)

write_csv(combined_regression_tbl, 'combined_regression_table.csv')


combined_regression_plot <- combined_regression_tbl |>
    ggplot(aes(y = energy_diff_peratom, x = qm9_u0)) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linetype = 'dashed') +
    geom_point(size = 0, alpha = 0.5) +
    geom_smooth(method = lmrob, se = FALSE) + 
    coord_obs_pred() +
    facet_wrap(~method, nrow = 1) +
    ylab('Energy Difference per Atom') + xlab('QM9 DFT') +
    # Limits based on the range of the QM9 data
    xlim(lower_limit - .2, upper_limit + .2) +
    ylim(lower_limit - .2, upper_limit + .2)
ggsave('dftbplus_tbmalt_comparison_facet.png', combined_regression_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Linear model relating the DFTB+ values to the tbmalt values
dftb_model <- lmrob(tbmalt_energy_diff ~ dftb_energy_diff, data = pre_regression_table)
