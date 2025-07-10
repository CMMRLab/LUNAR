# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 19:15:05 2025

@author: jdkem
"""


import numpy as np
import math


def generate_graphene(nx, ny, a=1.42):
    """
    Generate graphene sheet atom positions (x, y, z=0).
    
    Parameters:
        nx (int): Number of unit cells along x.
        ny (int): Number of unit cells along y.
        a (float): C–C bond length (default 1.42 Å).
    
    Returns:
        positions (list of tuple): List of (x, y, z=0) positions.
    """
    # Unit cell vectors
    a1 = np.array([np.sqrt(3)*a, 0])
    a2 = np.array([np.sqrt(3)*a/2, 3*a/2])
    
    # Basis atoms in unit cell
    basis = [
        np.array([0, 0]),
        np.array([0, a])
    ]
    
    positions = []
    for i in range(nx):
        for j in range(ny):
            # Lattice translation vector
            R = i * a1 + j * a2
            for b in basis:
                pos = R + b
                positions.append((pos[0], pos[1], 0))
                
    # Center positions
    positions = np.array(positions)
    center = positions.mean(axis=0)
    rel_positions = positions - center
    print(rel_positions)
    
    return rel_positions


def trim_hexagon1(positions, side_length):
    """
    Trim a graphene flake to a symmetrical hexagon.
    
    Parameters:
        positions (list of tuple): Atom positions (x, y).
        side_length (float): Side length of the hexagon.
    
    Returns:
        trimmed_positions (list of tuple): Positions within hexagon.
    """
    trimmed = []
    for x, y, z in positions:
        # Transform to axial (hex) coordinates
        q = (np.sqrt(3)/3 * x - 1/3 * y) / side_length
        r = (2/3 * y) / side_length
        s = -q - r
        if max(abs(q), abs(r), abs(s)) <= 1:
            trimmed.append((x, y, z))
    return np.array(trimmed)



def trim_hexagon2(positions, hex_radius):
    """
    Trim graphene atoms to form a symmetric hexagonal region.

    Parameters:
        positions (list of tuple): Atom positions (x, y, z).
        hex_radius (float): Hexagon distance from center to vertex.

    Returns:
        hex_positions (list of tuple): Atoms within the hexagonal region.
    """
    hex_positions = []
    for x, y, z in positions:
        # Compute axial hex coordinates
        q = (np.sqrt(3)/3 * x - 1/3 * y) / hex_radius
        r = (2/3 * y) / hex_radius
        s = -q - r
        if max(abs(q), abs(r), abs(s)) <= 1:
            hex_positions.append((x, y, z))
    return np.array(hex_positions)


def trim_hexagon3(positions, hex_radius):
    """
    Cut atoms to form a symmetric hexagonal graphene flake.

    Parameters:
        positions (list of tuple): Atom positions (x, y, z).
        hex_radius (float): Distance from center to hexagon vertex.

    Returns:
        flake_positions (list of tuple): Atoms within hexagonal flake.
    """
    flake_positions = []
    for x, y, z in positions:
        # Calculate the six edge directions (hexagonal condition)
        px, py = abs(x), abs(y)
        if (np.sqrt(3) * px + py <= np.sqrt(3) * hex_radius and
            np.sqrt(3) * px - py <= np.sqrt(3) * hex_radius and
            py <= hex_radius):
            flake_positions.append((x, y, z))
    return np.array(flake_positions)


def trim_hexagon4(positions, hex_radius):
    """
    Cut graphene atoms to form a symmetric hexagonal flake.

    Parameters:
        positions (list of tuple): Atom positions (x, y, z).
        hex_radius (float): Distance from center to hexagon vertex.

    Returns:
        hex_flake (list of tuple): Atoms within the hexagon.
    """
    hex_flake = []
    for x, y, z in positions:
        # Compute absolute coordinates
        px, py = abs(x), abs(y)
        # Symmetric hexagon condition
        if (py <= hex_radius and
            np.sqrt(3) * px + py <= np.sqrt(3) * hex_radius and
            np.sqrt(3) * px - py <= np.sqrt(3) * hex_radius):
            hex_flake.append((x, y, z))
    return np.array(hex_flake)

# Use law of cosines to compute the "zigzag span"
r0 = 1.42
n = 4
span = math.sqrt(r0**2 + r0**2 - 2*r0*r0*math.cos(math.radians(120)))
side_length = (n - 1.5)*span
print(n, side_length)


positions = generate_graphene(20, 20, a=1.42)
hex_positions = trim_hexagon1(positions, side_length)
hex_positions = trim_hexagon2(positions, side_length)
hex_positions = trim_hexagon3(positions, side_length)
hex_positions = trim_hexagon4(positions, side_length)





def write_lmp(filename, hex_positions):
    # Writing new file with bonds information
    with open(filename,'w') as f: 
        header = 'Testing graphene'
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
        f.write(f'{len(hex_positions)} atoms\n\n')
        f.write('1 atom types\n\n')
        
        buffer = 5.0
        xlo = np.min(hex_positions[:,0]) - buffer
        xhi = np.max(hex_positions[:,0]) + buffer
        ylo = np.min(hex_positions[:,1]) - buffer
        yhi = np.max(hex_positions[:,1]) + buffer
        zlo = np.min(hex_positions[:,2]) - buffer
        zhi = np.max(hex_positions[:,2]) + buffer
        f.write('{:<12.8f} {:<12.8f} xlo xhi\n'.format(xlo, xhi))
        f.write('{:<12.8f} {:<12.8f} ylo yhi\n'.format(ylo, yhi))
        f.write('{:<12.8f} {:<12.8f} zlo zhi\n'.format(zlo, zhi))
        
        f.write('\nMasses\n\n')
        f.write(' 1   12.01100  # C\n')
        
        f.write('\nAtoms # full\n\n')
        for n, (x, y, z) in enumerate(hex_positions, 1):
            molid = 1
            charge = 0.0
            atomtype = 1
            f.write('{:^6} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f}\n'.format(n, molid, atomtype, charge, x, y, z))
    return

#write_lmp('qdot.data', hex_positions)









import matplotlib.pyplot as plt
import numpy as np

def trim_hexagon(positions, side_length, center=(0.0, 0.0)):
    trimmed = []
    cx, cy = center
    for x, y in positions:
        # Shift to hexagon center
        x_shifted, y_shifted = x - cx, y - cy
        q = (np.sqrt(3)/3 * x_shifted - 1/3 * y_shifted) / side_length
        r = (2/3 * y_shifted) / side_length
        s = -q - r
        if max(abs(q), abs(r), abs(s)) <= 1:
            trimmed.append((x, y))
    return trimmed

# Example: generate some atom positions in a 2D grid
positions = [(x, y) for x in np.linspace(-25, 25, 50)
                      for y in np.linspace(-25, 25, 50)]

# Hexagon centered at (5, 5)
side_length = 10.0
center = (5.0, 0.0)
hex_flake = trim_hexagon(positions, side_length, center)

# Plot
x_all, y_all = zip(*positions)
x_hex, y_hex = zip(*hex_flake)

plt.figure(figsize=(6, 6))
plt.scatter(x_all, y_all, s=10, c='lightgray', label='All Positions')
plt.scatter(x_hex, y_hex, s=10, c='red', label='Trimmed Hexagon')
plt.plot(center[0], center[1], 'bo', label='Hexagon Center', markersize=8)
plt.legend()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Hexagon Trimmed from Grid')
plt.grid(True)
plt.show()

