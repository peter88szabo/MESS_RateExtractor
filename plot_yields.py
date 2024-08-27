import matplotlib.pyplot as plt
import numpy as np


def chemact_vs_thermal(single_plot, dataclass, file_path, delta_e_text, ax):
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
        ax.setp(line_2, markerfacecolor=line_2.get_color())

    # Customize the axis
    ax.set_xlabel('Temperature (K)', fontsize=14)
    ax.set_ylabel('Yield', fontsize=14)
    ax.set_xlim(250, 500)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle=':')
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Add the text annotations
    ax.text(310, 0.61, r"       Chemical Activation", fontsize=11, color='blue', bbox=dict(facecolor='white', edgecolor='white'))
    ax.text(310, 0.58, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ $\delta$-HPALD", fontsize=10, color='blue', bbox=dict(facecolor='white', edgecolor='white'))

    ax.text(310, 0.41, r"       Thermal Activation", fontsize=11, color='red', bbox=dict(facecolor='white', edgecolor='white'))
    ax.text(310, 0.38, r"Z,Z'-OH-allyl + O$_2$ $\rightarrow$ Peroxy-2 $\rightarrow$ $\delta$-HPALD", fontsize=10, color='red', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the energy level text
    ax.text(460, 0.30, delta_e_text, fontsize=12, color='black')
    ax.text(460, 0.70, r"Case-1", fontsize=14, color='black', bbox=dict(facecolor='white', edgecolor='white'))

    # Add the legend
    header_r_p1, = ax.plot([], [], color='none', label="")  # Invisible line for the header of R → P1
    header_r_w1, = ax.plot([], [], color='none', label="")  # Invisible line for the header of R → W1

    handles = [header_r_p1] + handles_r_p1 + [header_r_w1] + handles_r_w1
    labels = ["Chemical Act."] + labels_r_p1 + ["Thermal Act."] + labels_r_w1

    ax.legend(handles=handles, labels=labels, loc='center right', ncol=2)

    if single_plot:
        # Create the legend with two columns
        plt.legend(handles=handles, labels=labels, loc='center right', ncol=2)
        plt.savefig('yields_chemact_vs_thermal.jpeg', format='jpeg', dpi=600)
        plt.show()


def plot_combined_yields(dataclass, direc, files, delta_e_texts):
    # Create a figure with 3 subplots (stacked vertically)
    fig, axs = plt.subplots(3, 1, figsize=(12, 24))

    # Loop through each file and create a plot in a different subplot
    for i, (file, delta_e_text) in enumerate(zip(files, delta_e_texts)):
        file_path = direc + file
        single_plot=False
        chemact_vs_thermal(single_plot, dataclass, file_path, delta_e_text, axs[i])

    plt.tight_layout()
    plt.savefig('yields_chemact_vs_thermal_combined.jpeg', format='jpeg', dpi=600)
    plt.show()

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

    plt.text(460, 0.30, r"$\Delta$E = 200 cm$^{-1}$",fontsize=12, color='black')
    plt.text(460, 0.70, r"Case-1",fontsize=14, color='black', bbox=dict(facecolor='white', edgecolor='white'))


    # Add invisible plot objects to act as headers
    header_r_p1, = plt.plot([], [], color='none', label="")  # Invisible line for the header of R → P1
    header_r_w1, = plt.plot([], [], color='none', label="")  # Invisible line for the header of R → W1

    # Combine the headers and curve labels for the two columns
    handles = [header_r_p1] + handles_r_p1 + [header_r_w1] + handles_r_w1
    labels = ["Chemical Act."] + labels_r_p1 + ["Thermal Act."] + labels_r_w1

    # Create the legend with two columns
    plt.legend(handles=handles, labels=labels, loc='center right', ncol=2)
    plt.savefig('yields_chemact_vs_thermal.jpeg', format='jpeg', dpi=600)
    plt.show()
    #========================================================================================



