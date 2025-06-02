import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import KDTree

# Parameters
n_groups = 5
points_per_group = 30
linking_length = 3
colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta']

# Generate well-separated groups on a circle
angles = np.linspace(0, 2*np.pi, n_groups, endpoint=False)
radius = 12
centers = np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))

points = []
for center in centers:
    points.append(np.random.normal(loc=center, scale=2.0, size=(points_per_group, 2)))
points = np.vstack(points)

# FoF setup
tree = KDTree(points)
n_points = len(points)
visited = np.zeros(n_points, dtype=bool)
group_labels = -1 * np.ones(n_points, dtype=int)
group_id = 0

# Track step-by-step animation
frame_sequence = []

for i in range(n_points):
    if not visited[i]:
        queue = [i]
        while queue:
            idx = queue.pop(0)
            if visited[idx]:
                continue
            visited[idx] = True
            group_labels[idx] = group_id
            frame_sequence.append((idx, group_id))  # Log each step
            neighbors = tree.query_ball_point(points[idx], linking_length)
            for nbr in neighbors:
                if not visited[nbr]:
                    queue.append(nbr)
        group_id += 1

# Plot setup
fig, ax = plt.subplots()
sc = ax.scatter(points[:, 0], points[:, 1], c='lightgrey')
highlight, = ax.plot([], [], 'o', color='black', markersize=12, markerfacecolor='none', markeredgewidth=2)
ax.set_title("FoF Group Discovery Step-by-Step")
ax.set_aspect('equal')
ax.set_xlim(points[:, 0].min() - 5, points[:, 0].max() + 5)
ax.set_ylim(points[:, 1].min() - 5, points[:, 1].max() + 5)

# Color memory
final_colors = ['lightgrey'] * n_points

def animate(i):
    if i >= len(frame_sequence):
        return sc, highlight

    idx, gid = frame_sequence[i]
    color = colors[gid % len(colors)]
    final_colors[idx] = color

    # Update point colors
    sc.set_color(final_colors)

    # Highlight current search index
    highlight.set_data([points[idx, 0]], [points[idx, 1]])

    return sc, highlight

ani = animation.FuncAnimation(
    fig, animate, frames=len(frame_sequence) + 10, interval=100, blit=True, repeat=False
)

plt.show()
