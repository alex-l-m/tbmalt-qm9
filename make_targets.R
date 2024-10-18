library(tidyverse)
library(robustbase)

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
    # Total number of atoms
    mutate(n_atoms = n_H + n_C + n_N + n_O) |>
    # Join with DFTB results
    inner_join(dftb, by = 'mol_id') |>
    # Filter out NA rows
    filter(!is.na(tbmalt_energy), !is.na(qm9_u0), !is.na(dftb_energy)) |>
    # Subtract out atomic contributions with regression
    mutate(tbmalt_energy_ref = predict(lmrob(tbmalt_energy ~ n_H + n_C + n_N + n_O)),
           qm9_u0_ref = predict(lmrob(qm9_u0 ~ n_H + n_C + n_N + n_O)),
           tbmalt_energy_diff = tbmalt_energy - tbmalt_energy_ref,
           qm9_u0_diff = qm9_u0 - qm9_u0_ref)

# Reference the QM9 values to the tbmalt atomic references by adding the energy
# difference to the fitted values
tbmalt_target <- pre_regression_table |>
    mutate(tbmalt_target = tbmalt_energy_ref + qm9_u0_diff) |>
    select(mol_id, tbmalt_target)

write_csv(tbmalt_target, 'tbmalt_target.csv')
