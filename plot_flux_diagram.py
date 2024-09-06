def flux_diagram(dataclass, tempr, press):
    '''
    tempr must be a string tempr = "300"
    press must be a string tempr = "760"
    '''
    import plotly.graph_objects as go


    # Extract the temperature list
    temperatures = dataclass.temp

    # Find the index for 300 K
    if tempr in temperatures:
        index_300K = temperatures.index(tempr)
    else:
        print(f"{tempr} K temperature not found in the dataset.")
        exit()

    # Extract the yields for each branching reaction at 760 Torr and 300 K, including 'R'
    #Scale down the yields from the previous step
    from_R_to_W1 = dataclass.get_temp_dependent_yields("R", "W1", press)[index_300K]
    from_W1_to_W5 = dataclass.get_temp_dependent_yields("W1", "W5", press)[index_300K]
    from_W1_to_W2 = dataclass.get_temp_dependent_yields("W1", "W2", press)[index_300K]


    yields = {
        "R": {
            "W1": dataclass.get_temp_dependent_yields("R", "W1", press)[index_300K],
            "W2": 0.0,#dataclass.get_temp_dependent_yields("R", "W2", "760")[index_300K],
            "W5": 0.0,#dataclass.get_temp_dependent_yields("R", "W5", "760")[index_300K],
            "P1": dataclass.get_temp_dependent_yields("R", "P1", press)[index_300K],
            "P4": dataclass.get_temp_dependent_yields("R", "P4", press)[index_300K],
            "P9": dataclass.get_temp_dependent_yields("R", "P9", press)[index_300K]
        },
        "W1": {
            "W2": from_R_to_W1 * dataclass.get_temp_dependent_yields("W1", "W2", press)[index_300K],
            "W5": from_R_to_W1 * dataclass.get_temp_dependent_yields("W1", "W5", press)[index_300K],
            "P1": from_R_to_W1 * dataclass.get_temp_dependent_yields("W1", "P1", press)[index_300K],
            "P4": 0.0,#dataclass.get_temp_dependent_yields("W1", "P4", "760")[index_300K],
            "P9": 0.0,#dataclass.get_temp_dependent_yields("W1", "P9", "760")[index_300K]
        },
        "W2": {
            "W1": from_R_to_W1 * from_W1_to_W2 * dataclass.get_temp_dependent_yields("W2", "W1", press)[index_300K],
            "W5": from_R_to_W1 * 0.0, #dataclass.get_temp_dependent_yields("W2", "W5", "760")[index_300K],
            "P1": from_R_to_W1 * 0.0, #dataclass.get_temp_dependent_yields("W2", "P1", "760")[index_300K],
            "P4": from_R_to_W1 * from_W1_to_W2 * dataclass.get_temp_dependent_yields("W2", "P4", press)[index_300K],
            "P9": from_R_to_W1 * from_W1_to_W2 * dataclass.get_temp_dependent_yields("W2", "P9", press)[index_300K]
        },
        "W5": {
            "W1": from_R_to_W1 * from_W1_to_W5 * dataclass.get_temp_dependent_yields("W5", "W1", press)[index_300K],
            "W2": 0.0,# dataclass.get_temp_dependent_yields("W5", "W2", "760")[index_300K],
            "P1": from_R_to_W1 * from_W1_to_W5 * dataclass.get_temp_dependent_yields("W5", "P1", press)[index_300K],
            "P4": 0.0,#dataclass.get_temp_dependent_yields("W5", "P4", "760")[index_300K],
            "P9": 0.0#dataclass.get_temp_dependent_yields("W5", "P9", "760")[index_300K]
        },
    }

    # Define nodes (each unique species is a node, including R)
    #labels = ["R", "W1", "W2", "W5", "P1", "P4", "P9"]
    labels = ["Z,Z'-OH-allyl + O2", "Peroxy-2", "Peroxy-3", "Peroxy-10", "δ-HPALD", "δ-Acid", "β-HPALD"]

    # Define source and target indices based on the branching yields
    sources = [
        0, 0, 0, 0, 0, 0,  # R to W1, W2, W5, P1, P4, P9
        1, 1, 1, 1, 1,  # W1 to others
        2, 2, 2, 2, 2,  # W2 to others
        3, 3, 3, 3, 3,  # W5 to others
    ]

    targets = [
        1, 2, 3, 4, 5, 6,  # R to W1, W2, W5, P1, P4, P9
        2, 3, 4, 5, 6,  # W1 to others
        1, 3, 4, 5, 6,  # W2 to others
        1, 2, 4, 5, 6,  # W5 to others
    ]

    # The flow (values) corresponding to the yields
    values = [
        yields["R"]["W1"], yields["R"]["W2"], yields["R"]["W5"], yields["R"]["P1"], yields["R"]["P4"], yields["R"]["P9"],  # R to others
        yields["W1"]["W2"], yields["W1"]["W5"], yields["W1"]["P1"], yields["W1"]["P4"], yields["W1"]["P9"],  # W1 to others
        yields["W2"]["W1"], yields["W2"]["W5"], yields["W2"]["P1"], yields["W2"]["P4"], yields["W2"]["P9"],  # W2 to others
        yields["W5"]["W1"], yields["W5"]["W2"], yields["W5"]["P1"], yields["W5"]["P4"], yields["W5"]["P9"],  # W5 to others
    ]

    # Create a gradient of colors for the links to simulate directionality
    colors = [
        "rgba(31, 119, 180, 0.6)", "rgba(255, 127, 14, 0.6)", "rgba(44, 160, 44, 0.6)",
        "rgba(214, 39, 40, 0.6)", "rgba(148, 103, 189, 0.6)", "rgba(140, 86, 75, 0.6)"
    ]

    link_colors = [
        colors[i % len(colors)] for i in range(len(sources))
    ]

    # Build the Sankey diagram with gradient coloring
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=labels,
            color="blue",  # You can customize node colors if desired
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            arrowlen=25,  # Add arrowheads to the links
            color=link_colors  # Use gradient colors to emphasize directionality
        )
    )])

    fig.show()


