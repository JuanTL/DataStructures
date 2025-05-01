import matplotlib.pyplot as plt
import numpy as np

# Function to load rectangles from file
def load_rectangles(filename):
    rectangles = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                min_x, min_y, max_x, max_y = map(float, line.strip().split())
                rectangles.append((min_x, min_y, max_x, max_y))
    except FileNotFoundError:
        print(f"File {filename} not found. Skipping rectangles.")
    return rectangles

# Function to load points from file
def load_points(filename):
    points = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                x, y = map(float, line.strip().split())
                points.append((x, y))
    except FileNotFoundError:
        print(f"File {filename} not found. Skipping points.")
    return points

# Function to plot rectangles and points
def plot_spatial_data(rectangles, points):
    fig, ax = plt.subplots()
    
    # Plot rectangles
    for rect in rectangles:
        min_x, min_y, max_x, max_y = rect
        width = max_x - min_x
        height = max_y - min_y
        rect_patch = plt.Rectangle((min_x, min_y), width, height, linewidth=1, edgecolor='blue', facecolor='none')
        ax.add_patch(rect_patch)
    
    # Plot points
    if points:
        points_array = np.array(points)
        ax.scatter(points_array[:, 0], points_array[:, 1], color='red', s=50, label='Points')
    
    # Set plot properties
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Rectangles and Points')
    ax.legend()
    ax.grid(True)
    plt.axis('equal')  # Equal scaling on both axes
    plt.show()

# Main function
def main():
    # File paths
    rectangles_file = 'rectangles0.txt'
    points_file = 'points0.txt'
    
    # Load data
    rectangles = load_rectangles(rectangles_file)
    points = load_points(points_file)
    
    # Plot the data
    plot_spatial_data(rectangles, points)

if __name__ == '__main__':
    main()
