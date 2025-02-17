import numpy as np
import time
import argparse
import time

# Define the decorator to measure execution time
def timing_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Record the start time
        result = func(*args, **kwargs)  # Call the original function
        end_time = time.time()  # Record the end time
        print(f"Function '{func.__name__}' executed in {end_time - start_time:.6f} seconds.")
        return result  # Return the result of the function
    return wrapper

@timing_decorator
def generate_points_in_cube(cube_side, dimension, num_points, min_distance, max_iter=100000):
    """
    Generate random points in a cube with a minimum distance between points, 
    applying periodic boundary conditions (PBC) using the minimum image convention.

    Parameters:
    - cube_side (float): The length of the cube's side.
    - dimension (int): The number of spatial dimensions.
    - num_points (int): The number of points to generate.
    - min_distance (float): The minimum allowed distance between points.
    - max_iter (int): Maximum number of iterations.

    Returns:
    - A numpy array of points, or None if generation fails.
    """
    points_v = np.full((num_points, dimension), np.nan)  # Initialize points array
    accepted = 0  # Start without any accepted points
    iterations = 0  # Iteration counter

    while accepted < num_points:
        point = np.random.rand(dimension) * cube_side  # Generate a random point
        is_valid = True

        # Check the new point against all accepted points
        for i in range(accepted):
            # Apply minimum image convention to calculate the distance
            delta = point - points_v[i]
            delta -= cube_side * np.round(delta / cube_side)  # Minimum image convention
            distance = np.linalg.norm(delta)

            # Check if the distance violates the minimum distance constraint
            if distance < min_distance:
                is_valid = False
                break

        if is_valid:
            points_v[accepted] = point
            print(f"Point {accepted+1} of {num_points} generated and accepted ({iterations} iter of {max_iter} maximum iterations)")
            accepted += 1

        iterations += 1
        if iterations >= max_iter:
            print("Failed to generate points: maximum iterations exceeded!")
            return None

    return points_v

def save_to_xyz_with_composition(filename, points, composition, cell=None):
    """
    Save points to an .xyz file with specified composition.

    Parameters:
    - filename (str): Name of the .xyz file to save.
    - points (numpy array): Array of shape (n_points, dimension) containing point coordinates.
    - composition (dict): A dictionary where keys are element symbols (str) and values are their 
                          percentages in the composition (float). Percentages should sum to 1.0.
    """
    if not np.isclose(sum(composition.values()), 1.0):
        raise ValueError("The composition percentages must sum to 1.0.")

    with open(filename, "w") as file:
        num_points = len(points)
        file.write(f"{num_points}\n")  # Number of points
        # We use the comment line of the xyz to activate the pbc (pbc="T T T") and 
        # specify the format of the file (Properties=species:S:1:pos:R:3, i.e. first
        # column for element and 3 columns for coordinates)
        file.write('pbc="T T T" Properties=species:S:1:pos:R:3\n') 
        # If cell information is provided, add it as a comment line
        if cell is not None:
            # Ensure the cell is a tuple or list with 6 values
            if len(cell) != 6:
                raise ValueError("Cell should contain 6 values: [a, b, c, alpha, beta, gamma].")
            file.write(f"Cell: {' '.join(map(str, cell))}\n")

        # Create a list of elements based on their percentages
        elements = []
        for element, percentage in composition.items():
            count = int(num_points * percentage)
            elements.extend([element] * count)
        
        # Shuffle elements to distribute them randomly
        np.random.shuffle(elements)

        for point, element in zip(points, elements):
            coords = " ".join(f"{coord:.6f}" for coord in point)  # Format coordinates
            file.write(f"{element} {coords}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate random points in a cube.")
    parser.add_argument("--cube_side", type=float, default=10.0, help="Cube side length")
    parser.add_argument("--dimension", type=int, default=3, help="Number of spatial dimensions")
    parser.add_argument("--num_points", type=int, required=True, help="Number of points to generate")
    parser.add_argument("--min_distance", type=float, default=1.0, help="Minimum distance between points")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility")
    parser.add_argument("--max_iter", type=int, default=100000, help="Maximum iterations to attempt")
    parser.add_argument("--filename", type=str, default="output.xyz", help="File.xyz in which the coordinates for the atoms are saved")
    parser.add_argument("--list_of_elements", type=str, nargs='+', default=["Element"], help="List of elements")
    parser.add_argument("--list_of_compositions", type=float, nargs='+', default=[1.0], help="List of compositions for each element")
    
    args = parser.parse_args()

    # Check for possible errors
    if len(args.list_of_elements) != len(args.list_of_compositions):
        raise ValueError("Failed! Number of elements and number of compositions must be equal!")
    
    elem_and_composition = dict(zip(args.list_of_elements, args.list_of_compositions))
    
    # Set the random seed for reproducibility
    if args.seed is not None:
        print(f"Using provided seed: {args.seed}")
        np.random.seed(args.seed)
    else:
        args.seed = int(time.time())
        print(f"No seed provided. Using time-based seed: {args.seed}")
        np.random.seed(args.seed)

    # Loop through each argument and print the parameter name and value
    for param, value in args.__dict__.items():
        print(f"parameter {param} has been set to {value}")

    points = generate_points_in_cube(
        cube_side=args.cube_side,
        dimension=args.dimension,
        num_points=args.num_points,
        min_distance=args.min_distance,
        max_iter=args.max_iter,
    )

    if points is not None:
        save_to_xyz_with_composition(args.filename, points, elem_and_composition)
        print(f"Generated points printed in {args.filename}")
    else:
        print("Failed to generate the requested points within the maximum number of iterations.")

