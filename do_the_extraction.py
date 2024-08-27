from MESS_extractor     import ChemNetwork
from plot_yields        import chemact_vs_thermal, plot_combined_yields, singple_plot_chemact_vs_thermal 
from make_yield_tables  import generate_temp_dependent_yield_table
from plot_flux_diagram  import flux_diagram 

direc="Example/"
file="Case1.out"
file_path = direc + file
print(file_path)

ME = ChemNetwork(file_path)

# Table 1: R → W1, W2, W5, P1, P4, P9 at pressures 100, 760, 1400 torr
reaction_pairs_R_to_others = [("R", "W1"), ("R", "W2"), ("R", "W5"), ("R", "P1"), ("R", "P4"), ("R", "P9")]
pressures = ['100', '760', '1400']  # Pressures as strings
file_name_1 = "Table_R_to_others"
title_1 = "Yields for R to Other Molecules"

# Table 2: W1 → R, W2, W5 at pressures 100, 760, 1400 torr
reaction_pairs_W1_to_others = [("W1", "R"), ("W1", "W2"), ("W1", "W5"), ("W1", "P1")]
file_name_2 = "Table_W1_to_others"
title_2 = "Yields for W1 to Other Molecules"

# Table 3: W2 → W1, P4, P9 at pressures 100, 760, 1400 torr
reaction_pairs_W2_to_others = [("W2", "W1"), ("W2", "P4"), ("W2", "P9"), ("W2", "R")]
file_name_3 = "Table_W2_to_others"
title_3 = "Yields for W2 to Other Molecules"

# Generate and display the three tables
generate_temp_dependent_yield_table(ME, reaction_pairs_R_to_others, pressures, file_name_1, title_1)
generate_temp_dependent_yield_table(ME, reaction_pairs_W1_to_others, pressures, file_name_2, title_2)
generate_temp_dependent_yield_table(ME, reaction_pairs_W2_to_others, pressures, file_name_3, title_3)

     

direc = "/home/peter/Dropbox/Research_Leuven/MESS_kinetics/ZhenProject_2024/Case-1_from_ZZ-allyl+O2/After_Review_Extra_Analysis_2024_08_20/"
files = [
    "ExtraTemp_and_Pressure_Case1_starting_from_ZZ_allyl+O2_version-0.9_Truhlar_50-100_cutoff.out",  # File 1
    "dE300_ExtraTemp_and_Pressure_Case1_starting_from_ZZ_allyl+O2_version-0.9_Truhlar_50-100_cutoff.out",  # File 2
    "dE400_ExtraTemp_and_Pressure_Case1_starting_from_ZZ_allyl+O2_version-0.9_Truhlar_50-100_cutoff.out"  # File 3
]

delta_e_texts = [r"$\Delta$E = 200 cm$^{-1}$",
                r"$\Delta$E = 300 cm$^{-1}$",
                r"$\Delta$E = 400 cm$^{-1}$"]

#singple_plot_chemact_vs_thermal(ME, file_path)

plot_combined_yields(ME, direc, files, delta_e_texts)

#flux_diagram(ME, "300")
