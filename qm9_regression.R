library(tidyverse)
library(robustbase)
library(cowplot)
library(ggdark)
library(tidymodels)
library(glue)
theme_set(dark_mode(theme_cowplot(font_size = 24)))
#theme_set(theme_cowplot(font_size = 12))
library(patchwork)

dftb_results_infile <- 'ase_dftb_qm9.csv'

dftb <- read_csv(dftb_results_infile, col_types = 
cols(
  mol_id = col_character(),
  energy = col_double())) |>
    rename(dftb_energy = energy)

TO_EV = 27.211386246
pre_regression_table <- read_csv('regression_table.csv', col_types = cols(
    mol_id = col_character(),
    qm9_id = col_integer(),
    n_H = col_double(),
    n_C = col_double(),
    n_N = col_double(),
    n_O = col_double(),
    tbmalt_energy = col_double(),
    qm9_u0 = col_double()
)) |>
    # Convert qm9 energy to electronvolts
    mutate(qm9_u0 = qm9_u0 * TO_EV) |>
    # Also the tbmalt energy
    mutate(tbmalt_energy = tbmalt_energy * TO_EV) |>
    # Total number of atoms
    mutate(n_atoms = n_H + n_C + n_N + n_O) |>
    # Join with DFTB results
    inner_join(dftb, by = 'mol_id') |>
    # Filter out NA rows
    filter(!is.na(tbmalt_energy), !is.na(qm9_u0), !is.na(dftb_energy)) |>
    # Subtract out atomic contributions with regression
    mutate(tbmalt_energy_ref = predict(lmrob(tbmalt_energy ~ n_H + n_C + n_N + n_O)),
           qm9_u0_ref = predict(lmrob(qm9_u0 ~ n_H + n_C + n_N + n_O)),
           dftb_energy_ref = predict(lmrob(dftb_energy ~ n_H + n_C + n_N + n_O)),
           tbmalt_energy_diff = tbmalt_energy - tbmalt_energy_ref,
           qm9_u0_diff = qm9_u0 - qm9_u0_ref,
           dftb_energy_diff = dftb_energy - dftb_energy_ref)


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

dftb_regression_plot <- pre_regression_table |>
    ggplot(aes(y = dftb_energy_diff / n_atoms, x = qm9_u0_diff / n_atoms)) +
    # y=x line
    geom_abline(intercept = 0, slope = 1, color = 'blue', linetype = 'dashed') +
    geom_point(size = 0, alpha = 0.5) +
    geom_smooth(method = lmrob, se = FALSE) + 
    coord_obs_pred() +
    ylab('DFTB+') + xlab('QM9 DFT')

two_plots <- dftb_regression_plot + regression_plot + plot_layout(ncol = 2, nrow = 1)
ggsave('dftbplus_tbmalt_comparison.png', two_plots, width = unit(11.51, 'in'), height = unit(4.75, 'in'))
