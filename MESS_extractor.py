import re
import numpy as np

class ChemNetwork:
    def __init__(self, file_path):
        self.energy_unit = None
        self.pressure_unit = None
        self.pressure = []
        self.temp = []
        self.species = []
        self.barriers = []
        self.type = {} #Dictionary for the type of different stationary points on the PES
        self.energy = {} #Dictionary to save energetics of the different stationary points on the PES
        self.rate_high_press = {}  # Dictionary to store rates with double keys
        self.yield_high_press = {}  # Dictionary to store rates with double keys
        self.rate_press_depn = {}  # Dictionary to store rates with double keys
        self.yield_press_depn = {}  # Dictionary to store rates with double keys

        file_content = self._read_file(file_path)

        self._parse_energies(file_content)
        self._parse_high_pressure_rate_coefficients(file_content)
        self._get_pressure_list(file_content)
        self._parse_pressure_dependent_rate_coefficients(file_content)
        self._get_yields()
        self._get_high_pressure_yields()


    def _get_yields(self):
        """
        Calculate the yields (branching ratios) for each reaction at each pressure and temperature point.
        This function normalizes the reaction rates by the total rates for each starting species, while handling None values.
        """

        self.rate_press_depn_toall = {}  # To store the sum of rates for each species (to-all rate)
        self.yield_press_depn = {}  # To store the calculated yields

        for (react_from, react_to), rate_matrix in self.rate_press_depn.items():
            if react_from not in self.rate_press_depn_toall:
                self.rate_press_depn_toall[react_from] = [[0 for _ in range(len(rate_matrix[0]))] for _ in range(len(rate_matrix))]

            # Add the rates to the to-all rate for this species, ignoring None values
            for i in range(len(rate_matrix)):  # Loop over temperature
                for j in range(len(rate_matrix[0])):  # Loop over pressure
                    if rate_matrix[i][j] is not None:
                        self.rate_press_depn_toall[react_from][i][j] += rate_matrix[i][j]

        # Now calculate the yields (branching ratios)
        for (react_from, react_to), rate_matrix in self.rate_press_depn.items():
            # Initialize the yield matrix if it doesn't exist
            if (react_from, react_to) not in self.yield_press_depn:
                self.yield_press_depn[(react_from, react_to)] = [[0 for _ in range(len(rate_matrix[0]))] for _ in range(len(rate_matrix))]

            # Avoid division by zero by checking the denominator and handle None values in rate_matrix
            for i in range(len(rate_matrix)):  # Loop over temperature
                for j in range(len(rate_matrix[0])):  # Loop over pressure
                    if rate_matrix[i][j] is not None and self.rate_press_depn_toall[react_from][i][j] != 0:
                        self.yield_press_depn[(react_from, react_to)][i][j] = (
                            rate_matrix[i][j] / self.rate_press_depn_toall[react_from][i][j]
                        )
                    else:
                        self.yield_press_depn[(react_from, react_to)][i][j] = None  # Set yield to 0 if rate is None or total rate is zero

    def _get_high_pressure_yields(self):
        """
        Calculate the high-pressure yields (branching ratios) for each reaction at each temperature point.
        This function normalizes the high-pressure reaction rates by the total rates for each starting species.
        """
        self.rate_high_press_toall = {}  # To store the sum of high-pressure rates for each species (to-all rate)
        self.yield_high_press = {}  # To store the calculated high-pressure yields

        for (react_from, react_to), rate_array in self.rate_high_press.items():
            if react_from not in self.rate_high_press_toall:
                self.rate_high_press_toall[react_from] = [0 for _ in range(len(rate_array))]

            # Add the rates to the to-all rate for this species, ignoring None values
            for i in range(len(rate_array)):  # Loop over temperature
                if rate_array[i] is not None:
                    self.rate_high_press_toall[react_from][i] += rate_array[i]

        # Now calculate the high-pressure yields (branching ratios)
        for (react_from, react_to), rate_array in self.rate_high_press.items():
            # Initialize the yield array if it doesn't exist
            if (react_from, react_to) not in self.yield_high_press:
                self.yield_high_press[(react_from, react_to)] = [0 for _ in range(len(rate_array))]

            # Avoid division by zero by checking the denominator and handle None values in rate_array
            for i in range(len(rate_array)):  # Loop over temperature
                if rate_array[i] is not None and self.rate_high_press_toall[react_from][i] != 0:
                    self.yield_high_press[(react_from, react_to)][i] = (
                        rate_array[i] / self.rate_high_press_toall[react_from][i]
                    )
                else:
                    self.yield_high_press[(react_from, react_to)][i] = 0  # Set yield to 0 if rate is None or total rate is zero


    def _read_file(self, file_path):
        with open(file_path, 'r') as file:
            file_content = file.readlines()
        return file_content

    def _get_energy_unit(self, line):
        match = re.search(r'([a-zA-Z/\-\s]+)\s*\):$', line)
        if match:
            last_word = match.group(1)  # The last word before the closing parenthesis
        else:
            raise ValueError("Could not find the energy unit")
        return last_word

    def _get_pressure_list(self,lines):
        press_section_start = False
        press_from_here = False
        press_stop_here = False
        self.pressure = []

        for line in lines:
            if "Pressure-Species Rate Tables:" in line:
                press_section_start = True
            elif press_section_start and "P(torr)" in line:
                press_from_here = True
            elif press_section_start and press_from_here and not press_stop_here:
                if "O-O" in line:
                    # Stop processing pressure values after the first "O-O"
                    press_section_start = False
                    press_stop_here = True
                    break  # Stop the loop here

                data = line.split()
                pr = data[0]

                if pr.isdigit() and pr not in self.pressure:
                    self.pressure.append(pr)

    def _parse_energies(self, lines):
        # Extract energetic information from the first part of the file
        well_sec = False
        prod_sec = False
        well_to_bimol_sec = False
        well_to_well_sec = False
        for line in lines:
            if "Wells" in line:
                if self.energy_unit == None:
                    self.energy_unit = self._get_energy_unit(line)
                well_sec = True
                prod_sec = False
                well_to_bimol_sec = False
                well_to_well_sec = False
                what = 'well'
            elif "Bimolecular Products" in line:
                if self.energy_unit == None:
                    self.energy_unit = self._get_energy_unit(line)
                well_sec = False
                prod_sec = True
                well_to_bimol_sec = False
                well_to_well_sec = False
                what = 'product'
            elif "Well-to-Bimolecular Barriers" in line:
                if self.energy_unit == None:
                    self.energy_unit = self._get_energy_unit(line)
                well_sec = False
                prod_sec = False
                well_to_bimol_sec = True
                well_to_well_sec = False
                what = 'barrier_well_to_bimol'
            elif "Well-to-Well Barriers" in line:
                if self.energy_unit == None:
                    self.energy_unit = self._get_energy_unit(line)
                well_sec = False
                prod_sec = False
                well_to_bimol_sec = False
                well_to_well_sec = True
                what = 'barrier_well_to_well'
            elif "___________" in line:
                well_sec = False
                prod_sec = False
                well_to_bimol_sec = False
                well_to_well_sec = False
                what = 'None'
                break


            # Extract energy from wells and barriers
            if 'Name' not in line and (well_sec or prod_sec or well_to_bimol_sec or well_to_well_sec) and (any(char.isdigit() for char in line)):
                parts = line.split()
                self.energy[parts[0]] = float(parts[1])
                self.type[parts[0]] = what 
                if well_sec or prod_sec:
                    self.species.append(parts[0])
                if well_to_bimol_sec or well_to_well_sec:
                    self.barriers.append(parts[0])

    def _parse_high_pressure_rate_coefficients(self, lines):
        hp_rate_section_start = False
        reactions = []
        Temp_First = 0

        for line in lines:
            if "High Pressure Rate Coefficients (Temperature-Species Rate Tables)" in line:
                hp_rate_section_start = True
                continue
            elif "Capture/Escape Rate Coefficients" in line:
                hp_rate_section_start = False
                break

            if hp_rate_section_start and "T(K)" in line:
                Temp_First += 1
                reactions = line.split()[1:] #We discarding the "T(K)" and keep the species-to-species table
                reactions = [tuple(reaction.split('->')) for reaction in reactions] #making pairs reactions = [('W5', 'W1'), ('W5', 'W2'), ...]
                continue


            # Extract rate data and store using double keys (tuple as key)
            if hp_rate_section_start and any(char.isdigit() for char in line) and "->" not in line:
                data = line.split()
                temperature = data[0]
                rate_data = data[1:]

                if Temp_First == 1 and temperature not in self.temp:
                    self.temp.append(temperature)

                #from_species = rate_data[0]

                for i, rate_value in enumerate(rate_data):

                    react_from =  reactions[i][0]
                    react_to =  reactions[i][1]
                    #print(i, react_from, react_to, rate_value)

                    # Initialize the dictionary key if it doesn't exist
                    if (react_from, react_to) not in self.rate_high_press:
                        self.rate_high_press[(react_from, react_to)] = []  # Initialize with an empty list

                    if rate_value == '***':
                        self.rate_high_press[(react_from, react_to)].append(None)
                    else:
                        self.rate_high_press[(react_from, react_to)].append(float(rate_value))
                    continue

    def _filter_after_splitting_string(self, unfiltered_reactions):
        # Filter out pairs where the second element is empty or elements not in self.species
        reactions = []
        for pair in unfiltered_reactions:
            if len(pair) == 2 and pair[1]:  # Check if the pair has both elements and the second one is not empty
                if pair[0] in self.species and pair[1] in self.species:  # Ensure both elements are in self.species
                    reactions.append(pair)
        return reactions



    def _parse_pressure_dependent_rate_coefficients(self, lines):
        pd_rate_section_start = False
        reactions = []
        Temp_First = 0
        pres_unit_first = 0

        for line in lines:
            if "Temperature-Species Rate Tables:" in line:
                pd_rate_section_start = True
            elif pd_rate_section_start and "Pressure = " in line:
                current_pressure = line.split()[2]

                pres_unit_first += 1
                if pres_unit_first == 1:
                    self.pressure_unit = line.split()[3]
            elif pd_rate_section_start and "_______________________" in line:
                pd_rate_section_start = False
                break

            if pd_rate_section_start and "T(K)" in line:
                Temp_First += 1
                unfiltered_reactions = line.split()[1:] #We discarding the "T(K)" and keep the species-to-species table
                unfiltered_reactions = [tuple(reaction.split('->')) for reaction in unfiltered_reactions] #making pairs reactions = [('W5', 'W1'), ('W5', 'W2'), ...]
                reactions = self._filter_after_splitting_string(unfiltered_reactions)


            if 'Pressure' not in line and pd_rate_section_start and any(char.isdigit() for char in line) and "->" not in line:
                data = line.split()
                temperature = data[0]

                rate_data = data[1:(len(reactions) + 1)]

                if Temp_First == 1 and temperature not in self.temp:
                    self.temp.append(temperature)

                for i, rate_value in enumerate(rate_data):

                    # Extract rate data and store using double keys (tuple as key)
                    react_from =  reactions[i][0]
                    react_to =  reactions[i][1]

                    # Initialize the dictionary key if it doesn't exist
                    if (react_from, react_to) not in self.rate_press_depn:
                        self.rate_press_depn[(react_from, react_to)] = [[None for _ in range(len(self.pressure))] for _ in range(len(self.temp))]

                    rate_value = None if rate_value == '***' else float(rate_value)

                    temp_index = self.temp.index(temperature)
                    pressure_index = self.pressure.index(current_pressure)

                    self.rate_press_depn[(react_from, react_to)][temp_index][pressure_index] = rate_value


                continue

    def print_temp_dependent_rates(self, react_from, react_to, target_pressure):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.rate_press_depn:
            print(f"Error: Reaction pair {reaction_pair} not found in the data.")
            return

        react_mat = self.rate_press_depn[reaction_pair]

        if target_pressure not in self.pressure:
            print(f"Error: Pressure {target_pressure} torr not found in the list of pressures.")
            return

        pressure_index = self.pressure.index(target_pressure)

        print(f"\n{react_from} => {react_to} Rates at Pressure = {target_pressure} {self.pressure_unit}:")
        print(f"T(K)      Rate(s-1)")

        for temp_index, temp_row in enumerate(react_mat):
            rate = temp_row[pressure_index]  # Extract the rate for the given pressure
            print(f"{self.temp[temp_index]}       {rate}")


    def print_temp_dependent_yields(self, react_from, react_to, target_pressure):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.yield_press_depn:
            print(f"Error: Reaction pair {reaction_pair} not found in the data.")
            return

        react_mat = self.yield_press_depn[reaction_pair]

        if target_pressure not in self.pressure:
            print(f"Error: Pressure {target_pressure} torr not found in the list of pressures.")
            return

        pressure_index = self.pressure.index(target_pressure)

        print(f"\n{react_from} => {react_to} Yields at Pressure = {target_pressure} {self.pressure_unit}:")
        print(f"T(K)       Yield")

        for temp_index, temp_row in enumerate(react_mat):
            rate_yield = temp_row[pressure_index]  # Extract the rate for the given pressure
            print(f"{self.temp[temp_index]}       {rate_yield}")


    def get_temp_dependent_rates(self, react_from, react_to, target_pressure):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.rate_press_depn:
            raise ValueError(f"Reaction pair {reaction_pair} not found in the data.")

        if target_pressure not in self.pressure:
            print(f"Error: Pressure {target_pressure} torr not found in the list of pressures.")
            return

        pressure_index = self.pressure.index(target_pressure)

        react_mat = self.rate_press_depn[reaction_pair]

        rates = [react_mat[temp_index][pressure_index] for temp_index in range(len(self.temp))]

        return rates


    def get_temp_dependent_yields(self, react_from, react_to, target_pressure):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.yield_press_depn:
            raise ValueError(f"Reaction pair {reaction_pair} not found in the data.")

        if target_pressure not in self.pressure:
            print(f"Error: Pressure {target_pressure} torr not found in the list of pressures.")
            return

        pressure_index = self.pressure.index(target_pressure)

        react_mat = self.yield_press_depn[reaction_pair]

        rate_yield = [react_mat[temp_index][pressure_index] for temp_index in range(len(self.temp))]

        return rate_yield



    def print_pressure_dependent_rates(self, react_from, react_to, target_temperature):
        # Define the reaction pair
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.rate_press_depn:
            print(f"Error: Reaction pair {reaction_pair} not found in the data.")
            return

        if target_temperature not in self.temp:
            print(f"Error: Temperature {target_temperature} K not found in the list of temperatures.")
            return

        temp_index = self.temp.index(target_temperature)

        react_mat = self.rate_press_depn[reaction_pair]

        print(f"\n{react_from} => {react_to} Rates at Temperature = {target_temperature} K:")
        print(f"P(torr)   Rate(s-1)")

        for pressure_index, pressure_value in enumerate(self.pressure):
            rate = react_mat[temp_index][pressure_index]  # Extract the rate for the given temperature
            print(f"{pressure_value}       {rate}")

        if reaction_pair in self.rate_high_press:
            high_pressure_rate = self.rate_high_press[reaction_pair][temp_index]  # Extract the high-pressure rate
            print(f"O-O       {high_pressure_rate}")
        else:
            print(f"O-O       High-pressure rate data not available for {reaction_pair}")


    def print_pressure_dependent_yields(self, react_from, react_to, target_temperature):
        # Define the reaction pair
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.yield_press_depn:
            print(f"Error: Reaction pair {reaction_pair} not found in the data.")
            return

        if target_temperature not in self.temp:
            print(f"Error: Temperature {target_temperature} K not found in the list of temperatures.")
            return

        temp_index = self.temp.index(target_temperature)

        react_mat = self.yield_press_depn[reaction_pair]

        print(f"\n{react_from} => {react_to} Yields at Temperature = {target_temperature} K:")
        print(f"P(torr)     Yields")

        for pressure_index, pressure_value in enumerate(self.pressure):
            rate = react_mat[temp_index][pressure_index]  # Extract the rate for the given temperature
            print(f"{pressure_value}       {rate}")

        if reaction_pair in self.yield_high_press:
            high_pressure_yield = self.yield_high_press[reaction_pair][temp_index]  # Extract the high-pressure rate
            print(f"O-O       {high_pressure_yield}")
        else:
            print(f"O-O       High-pressure yield data not available for {reaction_pair}")

    def get_pressure_dependent_rates(self, react_from, react_to, target_temperature):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.rate_press_depn:
            raise ValueError(f"Reaction pair {reaction_pair} not found in the data.")

        if target_temperature not in self.temp:
            raise ValueError(f"Temperature {target_temperature} K not found in the list of temperatures.")

        temp_index = self.temp.index(target_temperature)

        react_mat = self.rate_press_depn[reaction_pair]

        rates = [react_mat[temp_index][pressure_index] for pressure_index in range(len(self.pressure))]

        return rates

    def get_pressure_dependent_yields(self, react_from, react_to, target_temperature):
        reaction_pair = (react_from, react_to)

        if reaction_pair not in self.yield_press_depn:
            raise ValueError(f"Reaction pair {reaction_pair} not found in the data.")

        if target_temperature not in self.temp:
            raise ValueError(f"Temperature {target_temperature} K not found in the list of temperatures.")

        temp_index = self.temp.index(target_temperature)

        react_mat = self.yield_press_depn[reaction_pair]

        rate_yield = [react_mat[temp_index][pressure_index] for pressure_index in range(len(self.pressure))]

        return rate_yield



if __name__ == "__main__":
    file_path = 'Example/Case1_ZZ_allyl+O2.out'

    ME = ChemNetwork(file_path)

    # Access the energetic value for W5
    W5 = ME.energy.get('W5')
    W1 = ME.energy.get('W1')
    W2 = ME.energy.get('W2')
    R  = ME.energy.get('R')
    P1 = ME.energy.get('P1')
    P4 = ME.energy.get('P4')
    P9 = ME.energy.get('P9')

    unit = ME.energy_unit

    print(f"\nSpecies: ")
    print(ME.species)
    print(f"\nBarriers: ")
    print(ME.barriers)
    print(f"\nTemperature: ")
    print(ME.temp)
    print(f"\nPressure: ")
    print(ME.pressure)


    print(f"\nEnergies are in", unit)
    print(f"E[W1]: {W1} ", ME.type.get('W1'))
    print(f"E[W2]: {W2} ", ME.type.get('W2'))
    print(f"E[W5]: {W5} ",  ME.type.get('W5'))
    print(f"E[R]:   {R}  ", ME.type.get('R'))
    print(f"E[P1]: {P1}  ", ME.type.get('P1'))
    print(f"E[P4]: {P4} ", ME.type.get('P4'))
    print(f"E[P9]: {P9}  ", ME.type.get('P9'))


#    print(f"\nHigh Pressure Rates:")
#    for (key1, key2), val in ME.rate_high_press.items():
#        print(key1, key2, val)


    print(f"\nW1-->W5 High-Pressure Rates:")
    for i in ME.rate_high_press[("W1", "W5")]:
        print(i)

    print(f"\nW1-->R High-Pressure Rates:")
    for i in ME.rate_high_press[("W1", "R")]:
        print(i)


    print(f"\nW5-->P1 High-Pressure Rates:")
    for i in ME.rate_high_press[("W5", "P1")]:
        print(i)

#---------------------------------------------------------

    ME.print_temp_dependent_rates("W1", "W5", "100")
    rates = ME.get_temp_dependent_rates("W1", "W5", "100")
    print("Rates for W1->W5 at 100 K:", rates)
    print("Yields")
    ME.print_temp_dependent_yields("W1", "W2", "100")
    ME.print_temp_dependent_yields("W1", "W5", "100")
    ME.print_temp_dependent_yields("W1", "R", "100")
    ME.print_temp_dependent_yields("W1", "P1", "100")
    ME.print_temp_dependent_yields("W1", "P4", "100")

    ME.print_temp_dependent_rates("W1", "R", "100")
    ME.print_temp_dependent_rates("W5", "P1", "200")
    ME.print_temp_dependent_rates("W5", "P1", "1400")

    ME.print_pressure_dependent_rates("W1", "W5", "250")
    rates = ME.get_pressure_dependent_rates("W1", "W5", "250")
    print("Rates for W1->W5 at 250 K:", rates)
    print()
    print("Yield")
    ME.print_pressure_dependent_yields("W1", "W5", "250")
    ME.print_pressure_dependent_yields("W5", "P1", "250")
    ME.print_pressure_dependent_yields("W2", "P4", "250")

    ME.print_pressure_dependent_rates("W1", "R", "300")
    ME.print_pressure_dependent_rates("W5", "P1", "300")
    ME.print_pressure_dependent_rates("W5", "P1", "300")
    ME.print_pressure_dependent_rates("W5", "P1", "580")
