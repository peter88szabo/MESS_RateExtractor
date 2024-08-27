import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def generate_temp_dependent_yield_table(dataclass, reaction_pairs, pressures, file_name, title):
    """
    Function to generate a temperature-dependent yield table and save it to Excel and display it as a figure.
    :param reaction_pairs: List of tuples representing reaction pairs (reactant, product)
    :param pressures: List of pressures at which to compute the yields (as strings)
    :param file_name: Name for the Excel file
    :param title: Title for the figure
    """
    # Get the temperature range from the dataclass class
    temperatures = dataclass.temp  # Assuming dataclass.temp contains the list of temperatures

    # Prepare the multi-level column headers
    columns = pd.MultiIndex.from_tuples(
        [(pressure, f"{reactant}->{product}") for pressure in pressures for reactant, product in reaction_pairs],
        names=["Pressure", "Reaction"]
    )

    # Initialize the DataFrame with NaN values
    df = pd.DataFrame(np.nan, index=temperatures, columns=columns)

    # Fill in the data for each temperature
    for pressure in pressures:
        for reactant, product in reaction_pairs:
            # Get the temperature-dependent yields for the current reactant->product at this pressure
            yields = dataclass.get_temp_dependent_yields(reactant, product, pressure)

            # Ensure yields is a list and matches the length of temperatures
            if isinstance(yields, list) and len(yields) == len(temperatures):
                for i, temp in enumerate(temperatures):
                    df.loc[temp, (pressure, f"{reactant}->{product}")] = yields[i]
            else:
                print(f"Error: Mismatch in yield data for {reactant}->{product} at {pressure} torr.")

    df.to_excel(f"{file_name}.xlsx")

