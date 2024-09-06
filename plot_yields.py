import matplotlib.pyplot as plt
import numpy as np


def W2_yields(dataclass, file_path, delta_e_text, ax):
    # Define the reactions we are interested in
    reaction_1 = ("W2", "R")
    reaction_2 = ("W2", "W1")
    reaction_3 = ("W2", "P4")
    reaction_4 = ("W2", "P9")

    temp_values = np.array([float(temp) for temp in dataclass.temp])

    handles_w2_p4 = []
    labels_w2_p4 = []

    handles_w2_p9 = []
    labels_w2_p9 = []

    handles_w2_r = []
    labels_w2_r = []

    handles_w2_w1 = []
    labels_w2_w1 = []

    #colors = ['y', 'g', 'k', 'r', 'b', 'c', 'm']
    colors = ['b', 'g', 'k', 'r', 'r', 'c', 'k']
    indices_to_loop = [0, 4, 6] #preselected pressures (100, 760, 1400) torr
    for pressure_index in indices_to_loop:
        # Get the pressure value at the current index
        pressure_value = dataclass.pressure[pressure_index]

        yields_for_reaction_1 = dataclass.yield_press_depn.get(reaction_1)
        yields_for_reaction_2 = dataclass.yield_press_depn.get(reaction_2)
        yields_for_reaction_3 = dataclass.yield_press_depn.get(reaction_3)
        yields_for_reaction_4 = dataclass.yield_press_depn.get(reaction_4)

        yield_values_1 = [yields_for_reaction_1[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_2 = [yields_for_reaction_2[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_3 = [yields_for_reaction_3[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_4 = [yields_for_reaction_4[temp_index][pressure_index] for temp_index in range(len(temp_values))]

        # Plot the yields for "W2 -> R" (open circles)
        line_1, = ax.plot(temp_values, yield_values_1, linestyle='--', marker='o', color=colors[pressure_index % len(colors)], markerfacecolor='none', label=f'{pressure_value} torr')
        handles_w2_r.append(line_1)
        labels_w2_r.append(f'{pressure_value} torr')

        # Plot the yields for "W2 -- > W1" (full triangles)
        line_2, = ax.plot(temp_values, yield_values_2, linestyle=':', marker='>', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_w2_w1.append(line_2)
        labels_w2_w1.append(f'{pressure_value} torr')
        plt.setp(line_2, markerfacecolor=line_2.get_color())

        # Plot the yields for "W1 --> P4 " (full circles)
        line_3, = ax.plot(temp_values, yield_values_3, linestyle='-', marker='o', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_w2_p4.append(line_3)
        labels_w2_p4.append(f'{pressure_value} torr')
        plt.setp(line_3, markerfacecolor=line_3.get_color())

        # Plot the yields for "W1 --> P9 " (full stars)
        line_4, = ax.plot(temp_values, yield_values_4, linestyle=':', marker='s', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_w2_p9.append(line_4)
        labels_w2_p9.append(f'{pressure_value} torr')
        plt.setp(line_4, markerfacecolor=line_4.get_color())


    # Customize the axis
    ax.set_xlabel('Temperature (K)', fontsize=17)
    ax.set_ylabel('Yield', fontsize=16)
    ax.set_xlim(250, 500)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle=':')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # Add the text annotations
    #ax.text(300, 0.63, r"       Chemical Activation", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.57, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))

    #ax.text(300, 0.43, r"       Thermal Activation", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.37, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the energy level text
    ax.text(450, 0.001, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    ax.text(420, 0.001, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(350, 0.001, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(320, 0.001, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the legend
    header_w2_p4, = ax.plot([], [], color='none', label="")  # Invisible line for the header
    header_w2_r, = ax.plot([], [], color='none', label="")  # Invisible line for the header
    header_w2_w1, = ax.plot([], [], color='none', label="")  # Invisible line for the header
    header_w2_p9, = ax.plot([], [], color='none', label="")  # Invisible line for the header

    handles = [header_w2_r] + handles_w2_r + [header_w2_w1] + handles_w2_w1 + [header_w2_p4] + handles_w2_p4 + [header_w2_p9] + handles_w2_p9
    labels = [r"3 $\rightarrow$ 1"] + labels_w2_r + [r"3 $\rightarrow$ 2"] + labels_w2_w1 + [r"3 $\rightarrow$ 5"] + labels_w2_p4 + [r"3 $\rightarrow$ 9"] + labels_w2_p9

    ax.legend(handles=handles, labels=labels, loc='lower right', ncol=4, fontsize=13)  # Set legend font size



def plot_combined_W2_yields(direc, files, delta_e_texts):
    from MESS_extractor     import ChemNetwork
    # Create a figure with 3 subplots (stacked vertically)
    #fig, axs = plt.subplots(3, 1, figsize=(12, 16))
    fig, axs = plt.subplots(len(files), 1, figsize=(12, len(files)*5.33))

    # Loop through each file and create a plot in a different subplot
    for i, (file, delta_e_text) in enumerate(zip(files, delta_e_texts)):
        file_path = direc + file
        dataclass = ChemNetwork(file_path)
        W2_yields(dataclass, file_path, delta_e_text, axs[i])
        # Set the y-axis to logarithmic scale
        axs[i].set_yscale('log')
        axs[i].set_ylim(1e-5, 1)


     #Only set the x-label for the last subplot
        if i < len(files) - 1:
            axs[i].set_xlabel('')  # Remove x-label for the upper subplots
        else:
            axs[i].set_xlabel('Temperature (K)', fontsize=16)

    # Add the text annotations above the top subplot, outside the figure area
    #fig.text(0.2, 0.97, r"Chemical Activation:", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    #fig.text(0.3, 0.98, r"2: Peroxy-2 ; 3: Peroxy-3; 10: Peroxy-10", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    fig.text(0.5, 0.97, r"1: Z,Z'-OH-allyl + O$_2$;   5: $\delta$-acid + OH;   9: $\beta$-HPALD + OH;   2: Peroxy-2;  3: Peroxy-3", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

#    fig.text(0.7, 0.97, r"Thermal Activation:", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
 #   fig.text(0.7, 0.95, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

    plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.05)  # Add space between the subplots
    plt.savefig('Case-2_W2_yields_combined_100_200_300cm1.jpeg', format='jpeg', dpi=600)
    #plt.savefig('Case-2_W2_yields_combined_ExcessEnergy_5kcal.jpeg', format='jpeg', dpi=600)
    #plt.show()




def W1_yields(dataclass, file_path, delta_e_text, ax):
    # Define the reactions we are interested in
    reaction_1 = ("W1", "R")
    reaction_2 = ("W1", "W2")
    reaction_3 = ("W1", "P1")
    reaction_4 = ("W1", "W5")

    temp_values = np.array([float(temp) for temp in dataclass.temp])

    handles_w1_p1 = []
    labels_w1_p1 = []

    handles_w1_r = []
    labels_w1_r = []

    handles_w1_w2 = []
    labels_w1_w2 = []

    handles_w1_w5 = []
    labels_w1_w5 = []


    #colors = ['y', 'g', 'k', 'r', 'b', 'c', 'm']
    colors = ['b', 'g', 'k', 'r', 'r', 'c', 'k']
    indices_to_loop = [0, 4, 6] #preselected pressures (100, 760, 1400) torr
    for pressure_index in indices_to_loop:
        # Get the pressure value at the current index
        pressure_value = dataclass.pressure[pressure_index]

        '''
        yields_for_reaction_1 = np.array(dataclass.yield_press_depn.get(reaction_1))
        yields_for_reaction_2 = np.array(dataclass.yield_press_depn.get(reaction_2))
        yields_for_reaction_3 = np.array(dataclass.yield_press_depn.get(reaction_3))
        yields_for_reaction_4 = np.array(dataclass.yield_press_depn.get(reaction_4))
        yields_for_reaction_3plus4 = yields_for_reaction_3 + yields_for_reaction_4 


        print(f"p = {pressure_value} torr and T = {temp_values[0]} K")
        print(f"yields of {reaction_1}: ", yields_for_reaction_1[0][4] )
        print(f"yields of {reaction_2}: ", yields_for_reaction_2[0][4])
        print(f"yields of {reaction_3}: ", yields_for_reaction_3[0][4])
        print(f"yields of {reaction_4}: ", yields_for_reaction_4[0][4])
        print(f"yields of {reaction_3}+{reaction_4}: ", yields_for_reaction_3plus4[0][4])

        yield_values_1 = [yields_for_reaction_1[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_2 = [yields_for_reaction_2[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_3 = [yields_for_reaction_3[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        '''

        # Retrieve the yield data and replace None with 0
        yields_for_reaction_1 = np.array([[0 if v is None else v for v in row] for row in dataclass.yield_press_depn.get(reaction_1)])
        yields_for_reaction_2 = np.array([[0 if v is None else v for v in row] for row in dataclass.yield_press_depn.get(reaction_2)])
        yields_for_reaction_3 = np.array([[0 if v is None else v for v in row] for row in dataclass.yield_press_depn.get(reaction_3)])
        yields_for_reaction_4 = np.array([[0 if v is None else v for v in row] for row in dataclass.yield_press_depn.get(reaction_4)])

        # Using indexing with NumPy arrays
        yield_values_1 = yields_for_reaction_1[:, pressure_index]
        yield_values_2 = yields_for_reaction_2[:, pressure_index]
        yield_values_3 = yields_for_reaction_3[:, pressure_index] + yields_for_reaction_4[:, pressure_index]

        print(f"p = {pressure_value} torr and T = {temp_values[7]} K")
        print(f"yields of {reaction_1}: ", yield_values_1[7])
        print(f"yields of {reaction_2}: ", yield_values_2[7])
        print(f"yields of {reaction_3}: ", yield_values_3[7])

        # Plot the yields for "W1 -> R" (open circles)
        line_1, = ax.plot(temp_values, yield_values_1, linestyle='--', marker='o', color=colors[pressure_index % len(colors)], markerfacecolor='none', label=f'{pressure_value} torr')
        handles_w1_r.append(line_1)
        labels_w1_r.append(f'{pressure_value} torr')

        # Plot the yields for "W1 -- > W2" (full triangles)
        line_2, = ax.plot(temp_values, yield_values_2, linestyle=':', marker='>', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_w1_w2.append(line_2)
        labels_w1_w2.append(f'{pressure_value} torr')
        plt.setp(line_2, markerfacecolor=line_2.get_color())

        # Plot the yields for "W1 --> W5 + P1 together" (full circles)
        line_3, = ax.plot(temp_values, yield_values_3, linestyle='-', marker='o', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_w1_p1.append(line_3)
        labels_w1_p1.append(f'{pressure_value} torr')
        plt.setp(line_3, markerfacecolor=line_3.get_color())

    # Customize the axis
    ax.set_xlabel('Temperature (K)', fontsize=17)
    ax.set_ylabel('Yield', fontsize=16)
    ax.set_xlim(250, 450)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle=':')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # Add the text annotations
    #ax.text(300, 0.63, r"       Chemical Activation", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.57, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))

    #ax.text(300, 0.43, r"       Thermal Activation", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.37, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the energy level text
    ax.text(400, 0.30, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    ax.text(370, 0.30, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(260, 0.01, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(260, 0.02, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the legend
    header_w1_p1, = ax.plot([], [], color='none', label="")  # Invisible line for the header
    header_w1_r, = ax.plot([], [], color='none', label="")  # Invisible line for the header
    header_w1_w2, = ax.plot([], [], color='none', label="")  # Invisible line for the header

    handles = [header_w1_r] + handles_w1_r + [header_w1_p1] + handles_w1_p1 + [header_w1_w2] + handles_w1_w2
    labels = [r"2 $\rightarrow$ 1"] + labels_w1_r + [r"2 $\rightarrow$ 10+11"] + labels_w1_p1 + [r"2 $\rightarrow$ 3"] + labels_w1_w2

    ax.legend(handles=handles, labels=labels, loc='upper left', ncol=3, fontsize=13)  # Set legend font size



def plot_combined_W1_yields(direc, files, delta_e_texts):
    from MESS_extractor     import ChemNetwork
    # Create a figure with 3 subplots (stacked vertically)
    #fig, axs = plt.subplots(3, 1, figsize=(12, 16))
    fig, axs = plt.subplots(len(files), 1, figsize=(12, len(files)*5.33))

    # Loop through each file and create a plot in a different subplot
    for i, (file, delta_e_text) in enumerate(zip(files, delta_e_texts)):
        file_path = direc + file
        dataclass = ChemNetwork(file_path)
        W1_yields(dataclass, file_path, delta_e_text, axs[i])
        # Set the y-axis to logarithmic scale
        axs[i].set_yscale('log')
        axs[i].set_ylim(1e-4, 1)


     #Only set the x-label for the last subplot
        if i < len(files) - 1:
            axs[i].set_xlabel('')  # Remove x-label for the upper subplots
        else:
            axs[i].set_xlabel('Temperature (K)', fontsize=16)

    # Add the text annotations above the top subplot, outside the figure area
    #fig.text(0.2, 0.97, r"Chemical Activation:", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    #fig.text(0.3, 0.98, r"2: Peroxy-2 ; 3: Peroxy-3; 10: Peroxy-10", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    fig.text(0.5, 0.97, r"1: Z,Z'-OH-allyl + O$_2$;   11: $\delta$-HPALD + OH;   2: Peroxy-2;  3: Peroxy-3;  10: Peroxy-10", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

#    fig.text(0.7, 0.97, r"Thermal Activation:", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
 #   fig.text(0.7, 0.95, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

    plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.05)  # Add space between the subplots
    plt.savefig('Case-2_W1_yields_combined_100_200_300cm1.jpeg', format='jpeg', dpi=600)
    #plt.savefig('Case-2_W1_yields_combined_ExcessEnergy_5kcal.jpeg', format='jpeg', dpi=600)
    #plt.show()


def chemact_vs_thermal(dataclass, file_path, delta_e_text, ax):
    # Define the reactions we are interested in
    reaction_1 = ("R", "W1")
    reaction_2 = ("R", "P1")

    temp_values = np.array([float(temp) for temp in dataclass.temp])

    handles_r_p1 = []
    labels_r_p1 = []
    handles_r_w1 = []
    labels_r_w1 = []

    colors = ['y', 'g', 'k', 'r', 'b', 'c', 'm']

    # Loop through each pressure and plot the yield for both reactions
    for pressure_index, pressure_value in enumerate(dataclass.pressure):
        # Extract yields for both reactions at the current pressure
        yields_for_reaction_1 = dataclass.yield_press_depn.get(reaction_1)
        yields_for_reaction_2 = dataclass.yield_press_depn.get(reaction_2)

        yield_values_1 = [yields_for_reaction_1[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_2 = [yields_for_reaction_2[temp_index][pressure_index] for temp_index in range(len(temp_values))]

        # Plot the yields for "R --> W1" (open circles)
        line_1, = ax.plot(temp_values, yield_values_1, linestyle='--', marker='o', color=colors[pressure_index % len(colors)], markerfacecolor='none', label=f'{pressure_value} torr')
        handles_r_w1.append(line_1)
        labels_r_w1.append(f'{pressure_value} torr')

        # Plot the yields for "R --> P1" (full circles)
        line_2, = ax.plot(temp_values, yield_values_2, linestyle='-', marker='o', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_r_p1.append(line_2)
        labels_r_p1.append(f'{pressure_value} torr')
        plt.setp(line_2, markerfacecolor=line_2.get_color())

    # Customize the axis
    ax.set_xlabel('Temperature (K)', fontsize=17)
    ax.set_ylabel('Yield', fontsize=16)
    ax.set_xlim(250, 500)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle=':')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # Add the text annotations
    #ax.text(300, 0.63, r"       Chemical Activation", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.57, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD", fontsize=13, color='blue', bbox=dict(facecolor='white', edgecolor='white'))

    #ax.text(300, 0.43, r"       Thermal Activation", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(300, 0.37, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=13, color='red', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the energy level text
    ax.text(450, 0.80, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    ax.text(420, 0.80, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(350, 0.80, delta_e_text, fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))
    #ax.text(320, 0.80, r"Case-2", fontsize=17, color='black', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the legend
    header_r_p1, = ax.plot([], [], color='none', label="")  # Invisible line for the header of R → P1
    header_r_w1, = ax.plot([], [], color='none', label="")  # Invisible line for the header of R → W1

    handles = [header_r_p1] + handles_r_p1 + [header_r_w1] + handles_r_w1
    labels = ["Chemical Act."] + labels_r_p1 + ["Thermal Act."] + labels_r_w1

    ax.legend(handles=handles, labels=labels, loc='center right', ncol=2, fontsize=13)  # Set legend font size



def plot_combined_yields(direc, files, delta_e_texts):
    from MESS_extractor     import ChemNetwork
    # Create a figure with 3 subplots (stacked vertically)
    #fig, axs = plt.subplots(3, 1, figsize=(12, 16))
    fig, axs = plt.subplots(len(files), 1, figsize=(12, len(files)*5.33))

    # Loop through each file and create a plot in a different subplot
    for i, (file, delta_e_text) in enumerate(zip(files, delta_e_texts)):
        file_path = direc + file
        dataclass = ChemNetwork(file_path)
        chemact_vs_thermal(dataclass, file_path, delta_e_text, axs[i])

     #Only set the x-label for the last subplot
        if i < len(files) - 1:
            axs[i].set_xlabel('')  # Remove x-label for the upper subplots
        else:
            axs[i].set_xlabel('Temperature (K)', fontsize=16)

    # Add the text annotations above the top subplot, outside the figure area
    fig.text(0.3, 0.97, r"Chemical Activation:", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    fig.text(0.3, 0.95, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD", fontsize=15, color='blue', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

    fig.text(0.7, 0.97, r"Thermal Activation:", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
    fig.text(0.7, 0.95, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=15, color='red', ha='center', bbox=dict(facecolor='white', edgecolor='white'))

    plt.subplots_adjust(hspace=0.15, top=0.94, bottom=0.05)  # Add space between the subplots
    plt.savefig('Case-2_yields_chemact_vs_thermal_combined_100_200_300cm1.jpeg', format='jpeg', dpi=600)
    #plt.savefig('Case-2_yields_chemact_vs_thermal_combined_ExcessEnergy_5kcal.jpeg', format='jpeg', dpi=600)
    #plt.show()

def singple_plot_chemact_vs_thermal(dataclass, file_path):

    # Define the reactions we are interested in
    reaction_1 = ("R", "W1")
    reaction_2 = ("R", "P1")

    temp_values = np.array([float(temp) for temp in dataclass.temp])

    plt.figure(figsize=(12, 8))

    handles = []
    labels = []
    handles_r_p1 = []
    labels_r_p1 = []

    handles_r_w1 = []
    labels_r_w1 = []

    #colors = ['y', 'g', 'm', 'c', 'b', 'r', 'k']
    colors = ['y', 'g', 'k', 'r', 'b', 'c', 'm']
    # Loop through each pressure and plot the yield for both reactions
    for pressure_index, pressure_value in enumerate(dataclass.pressure):
        # Extract yields for both reactions at the current pressure
        yields_for_reaction_1 = dataclass.yield_press_depn.get(reaction_1)
        yields_for_reaction_2 = dataclass.yield_press_depn.get(reaction_2)

        yield_values_1 = [yields_for_reaction_1[temp_index][pressure_index] for temp_index in range(len(temp_values))]
        yield_values_2 = [yields_for_reaction_2[temp_index][pressure_index] for temp_index in range(len(temp_values))]


        # Plot the yields for "R --> W1" (open circles) with a unique label for the legend
        line_1, = plt.plot(temp_values, yield_values_1, linestyle='--', marker='o', color=colors[pressure_index % len(colors)], markerfacecolor='none', label=f'{pressure_value} torr')
        handles_r_w1.append(line_1)
        labels_r_w1.append(f'{pressure_value} torr')
        # Make the marker face color match the line color

        # Plot the yields for "R --> P1" (full circles) with a unique label for the legend
        line_2, = plt.plot(temp_values, yield_values_2, linestyle='-', marker='o', color=colors[pressure_index % len(colors)], label=f'{pressure_value} torr')
        handles_r_p1.append(line_2)
        labels_r_p1.append(f'{pressure_value} torr')
        # Make the marker face color match the line color
        plt.setp(line_2, markerfacecolor=line_2.get_color())

    #plt.title(f'Yields for Reactions R --> W1 and R --> P1 as a Function of Temperature')
    plt.xlabel('Temperature (K)',fontsize=14)
    plt.ylabel('Yield',fontsize=14)
    plt.xlim(250,500)
    plt.ylim(0,1)
    plt.grid(True,linestyle=':')
    plt.rcParams.update({'font.size': 12}) 
    plt.xticks(fontsize=14)  # X-axis numbers
    plt.yticks(fontsize=14)

    plt.text(310, 0.61, r"       Chemical Activation",fontsize=11, color='blue', bbox=dict(facecolor='white', edgecolor='white'))
    plt.text(310, 0.58, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD",fontsize=10, color='blue', bbox=dict(facecolor='white', edgecolor='white'))

    plt.text(310, 0.41, r"       Thermal Activation",fontsize=11, color='red', bbox=dict(facecolor='white', edgecolor='white'))
    plt.text(310, 0.38, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD",fontsize=10, color='red', bbox=dict(facecolor='white', edgecolor='white'))

    plt.text(450, 0.70, r"$\Delta$E = 200 cm$^{-1}$",fontsize=12, color='black')
    plt.text(420, 0.70, r"Case-2",fontsize=14, color='black', bbox=dict(facecolor='white', edgecolor='white'))


    # Add invisible plot objects to act as headers
    header_r_p1, = plt.plot([], [], color='none', label="")  # Invisible line for the header of R → P1
    header_r_w1, = plt.plot([], [], color='none', label="")  # Invisible line for the header of R → W1

    # Combine the headers and curve labels for the two columns
    handles = [header_r_p1] + handles_r_p1 + [header_r_w1] + handles_r_w1
    labels = ["Chemical Act."] + labels_r_p1 + ["Thermal Act."] + labels_r_w1

    # Create the legend with two columns
    plt.legend(handles=handles, labels=labels, loc='center right', ncol=2)
    plt.savefig('Case-2_yields_chemact_vs_thermal.jpeg', format='jpeg', dpi=600)
    plt.show()
    #========================================================================================



