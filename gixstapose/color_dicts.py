import fresnel


cpk_colors = {
    "H": fresnel.color.linear([1.00, 1.00, 1.00]),  # white
    "C": fresnel.color.linear([0.30, 0.30, 0.30]),  # grey
    "N": fresnel.color.linear([0.13, 0.20, 1.00]),  # dark blue
    "O": fresnel.color.linear([1.00, 0.13, 0.00]),  # red
    "F": fresnel.color.linear([0.12, 0.94, 0.12]),  # green
    "Cl": fresnel.color.linear([0.12, 0.94, 0.12]),  # green
    "Br": fresnel.color.linear([0.60, 0.13, 0.00]),  # dark red
    "I": fresnel.color.linear([0.40, 0.00, 0.73]),  # dark violet
    "He": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Ne": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Ar": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Xe": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Kr": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "P": fresnel.color.linear([1.00, 0.60, 0.00]),  # orange
    "S": fresnel.color.linear([1.00, 0.90, 0.13]),  # yellow
    "B": fresnel.color.linear([1.00, 0.67, 0.47]),  # peach
    "Li": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Na": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "K": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Rb": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Cs": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Fr": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Be": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Mg": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ca": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Sr": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ba": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ra": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ti": fresnel.color.linear([0.60, 0.60, 0.60]),  # grey
    "Fe": fresnel.color.linear([0.87, 0.47, 0.00]),  # dark orange
    "default": fresnel.color.linear([0.87, 0.47, 1.00]),  # pink
}

bsu_colors = [
        fresnel.color.linear([0.00, 0.20, 0.63]), # blue
        fresnel.color.linear([0.84, 0.30, 0.04]), # orange
        fresnel.color.linear([0.00, 0.23, 0.44]), # blue
        fresnel.color.linear([0.00, 0.42, 0.65]), # blue
        fresnel.color.linear([0.00, 0.45, 0.81]), # blue
        fresnel.color.linear([0.25, 0.38, 0.60]), # blue
        fresnel.color.linear([0.17, 0.99, 1.00]), # blue
        fresnel.color.linear([1.00, 0.39, 0.03]), # orange
        fresnel.color.linear([1.00, 0.40, 0.00]), # orange
        ]

# Made space to add more later
radii_dict = {"H": 0.05, "default": 0.06}
