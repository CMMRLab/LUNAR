This examples folder holds some LAMMPS .log.lammps files for testing the log_analysis.py GUI. The
example .log.lammps files are from a vareity of Josh's projects where the exact data is not 
important, but just serves to help users of log_analysis.py with some quick to load files.

This directory also contains a "using read_log.py for personal use" directory with the following
two python scripts to show how to use read_log.py for your own .log.lammps analysis in case the
log_analysis.py GUI can not provide enough analysis:
  - read_logfile_and_plot.py  which shows how to read data, load data, and plot data for .log.lammps files
  - read_logfile_and_write_2_csv which shows how to read data, load data, and write data to .csv for plotting in excel








INFO FOR JOSH TO KEEP TRACK OF WHERE THE EXAMPLE FILES CAME FROM (OTHERS PLEASE IGNORE THIS SECTION):
relax_rep_1_pxld_60_IFF.log.lammps -> property=density_ts=0.5.log.lammps

tensile_L10_rep_1_1_full_xld_erate_2e8_IFF-RX0.log.lammps -> property=tensile_modulus_x_strain_rate=2e8.log.lammps
tensile_L10_rep_1_2_full_xld_erate_2e8_IFF-RX0.log.lammps -> property=tensile_modulus_y_strain_rate=2e8.log.lammps
tensile_L10_rep_1_3_full_xld_erate_2e8_IFF-RX0.log.lammps -> property=tensile_modulus_z_strain_rate=2e8.log.lammps

shear_L10_rep_2_1_full_xld_erate_2e8_IFF-RX0_mmorse.log.lammps -> property=shear_modulus_xy_strain_rate=2e8.log.lammps
shear_L10_rep_1_2_full_xld_erate_2e8_IFF-RX0_mmorse.log.lammps -> property=shear_modulus_xz_strain_rate=2e8.log.lammps
shear_L10_rep_1_3_full_xld_erate_2e8_IFF-RX0_mmorse.log.lammps -> property=shear_modulus_yz_strain_rate=2e8.log.lammps

BulkM_EPON_862_rep_1_IFF.log.lammps -> property=bulk_modulus_v0_v1_calc_ts=0.5_p0=1atm_p1=5000atm.log.lammps

Tg_CTE_rep_1_mech_hrate_50_temps_100-770.log.lammps -> properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps