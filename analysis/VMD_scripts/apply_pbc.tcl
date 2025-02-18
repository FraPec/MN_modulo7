# Set up the periodic boundary conditions
pbc set {14.8 14.8 14.8 90 90 90} -all  ;# Define the PBC dimensions
pbc box                                      ;# Initialize the box for visualization
pbc wrap -all                                ;# Wrap atoms into the primary box
