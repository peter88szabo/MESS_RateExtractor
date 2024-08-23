import re

class ChemNetwork:
    def __init__(self, file_path):
        self.energy_unit = None
        self.temp = []
        self.species = []
        self.barriers = []
        self.type = {} #Dictionary for the type of different stationary points on the PES
        self.energy = {} #Dictionary to save energetics of the different stationary points on the PES
        self.rate_high_press = {}  # Dictionary to store rates with double keys
        self.rate_press_depn = {}  # Dictionary to store rates with double keys

        file_content = self._read_file(file_path)

        self._parse_energies(file_content)
        self._parse_high_pressure_rate_coefficients(file_content)
        #self._parse_pressure_dependent_rate_coefficients(file_content)

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

                if Temp_First == 1:
                    self.temp.append(temperature)

                #from_species = rate_data[0]

                for i, rate_value in enumerate(rate_data):

                    react_from =  reactions[i][0]
                    react_to =  reactions[i][1]
                    print(i, react_from, react_to, rate_value)

                    # Initialize the dictionary key if it doesn't exist
                    if (react_from, react_to) not in self.rate_high_press:
                        self.rate_high_press[(react_from, react_to)] = []  # Initialize with an empty list

                    if rate_value == '***':
                        self.rate_high_press[(react_from, react_to)].append(None)
                    else:
                        self.rate_high_press[(react_from, react_to)].append(float(rate_value))
                    continue




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
    for spec in ME.species:
        print(spec)
    print(f"\nBarriers: ")
    for barr in ME.barriers:
        print(barr)
    print(f"\nTemperature: ")
    for temp in ME.temp:
        print(temp)

    print(f"\nEnergies are in", unit)
    print(f"E[W1]: {W1} ", ME.type.get('W1'))
    print(f"E[W2]: {W2} ", ME.type.get('W2'))
    print(f"E[W5]: {W5} ",  ME.type.get('W5'))
    print(f"E[R]:   {R}  ", ME.type.get('R'))
    print(f"E[P1]: {P1}  ", ME.type.get('P1'))
    print(f"E[P4]: {P4} ", ME.type.get('P4'))
    print(f"E[P9]: {P9}  ", ME.type.get('P9'))


    print(f"\nHigh Pressure Rates:")
    for (key1, key2), val in ME.rate_high_press.items():
        print(key1, key2, val)


    print(f"\nW1-->W5 Rates:")
    for i in ME.rate_high_press[("W1", "W5")]:
        print(i)

    print(f"\nW1-->R Rates:")
    for i in ME.rate_high_press[("W1", "R")]:
        print(i)


    print(f"\nW5-->P1 Rates:")
    for i in ME.rate_high_press[("W5", "P1")]:
        print(i)

    # Access the high-pressure rate for W1 -> W5
    #rate = ME.rate_high_press.get(('W1', 'W5'))
    #if rate is not None:
    #    print(f"The high-pressure rate for 'W1' -> 'W5' is: {rate}")
    #else:
    #    print("No high-pressure rate found for 'W1' -> 'W5'")

