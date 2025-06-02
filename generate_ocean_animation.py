#!/usr/bin/env python3

import subprocess
import os
import time
import numpy as np
import imageio
import tempfile
import trimesh
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

def run_ocean_generator(time_val, output_path, width=20, length=20, width_segments=128, 
                       length_segments=128, amplitude=1.5, choppiness=0.5, wind_speed=5.0):
    """Run the Rust ocean generator with the specified parameters"""
    cmd = [
        "./target/debug/ocean-generator",
        "-w", str(width),
        "-l", str(length),
        "--width-segments", str(width_segments),
        "--length-segments", str(length_segments),
        "--wind-speed", str(wind_speed),
        "-a", str(amplitude),
        "-c", str(choppiness),
        "-t", str(time_val),
        "-o", output_path
    ]
    
    print(f"Generating ocean surface at time={time_val}...")
    subprocess.run(cmd, check=True)
    return output_path

def render_glb_to_image(glb_path, output_image_path, width=800, height=600, dpi=100):
    """Render a GLB file to an image using matplotlib"""
    try:
        # Load the GLB with trimesh - it returns a Scene object
        scene = trimesh.load(glb_path)
        
        # Get the mesh from the scene
        mesh = None
        if isinstance(scene, trimesh.Scene):
            # If we have a scene, extract the first mesh
            for geometry_name in scene.geometry:
                current_mesh = scene.geometry[geometry_name]
                if isinstance(current_mesh, trimesh.Trimesh):
                    mesh = current_mesh
                    break
        else:
            # If trimesh.load returned a mesh directly
            mesh = scene
            
        if mesh is None or not isinstance(mesh, trimesh.Trimesh):
            print(f"Error: Could not find a valid mesh in {glb_path}")
            # Create a simple placeholder image
            img = Image.new('RGB', (width, height), color=(100, 149, 237))  # Cornflower blue
            img.save(output_image_path)
            return output_image_path
            
        # Create a figure and 3D axis
        fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
        ax = fig.add_subplot(111, projection='3d')
        
        # Get vertices and faces
        vertices = mesh.vertices
        faces = mesh.faces
        
        # Get the z values for coloring
        z = vertices[:, 2]
        
        # Plot the 3D surface using triangulation
        surf = ax.plot_trisurf(
            vertices[:, 0], vertices[:, 1], vertices[:, 2],
            triangles=faces,
            cmap=cm.viridis,  # Using viridis for better color contrast
            linewidth=0,
            antialiased=True,
            alpha=0.9
        )
        
        # Add a color bar which maps values to colors
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Height (Z)')
        
        # Set axis labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Set viewing angle for better visualization of waves
        ax.view_init(30, 45)  # elevation, azimuth
        
        # Ensure the aspect ratio is equal for all axes
        # Exaggerate the z-axis to make waves more visible
        max_range = max([vertices[:, 0].max() - vertices[:, 0].min(),
                         vertices[:, 1].max() - vertices[:, 1].min(),
                         vertices[:, 2].max() - vertices[:, 2].min() * 2])  # Amplify z range
        
        mid_x = (vertices[:, 0].max() + vertices[:, 0].min()) * 0.5
        mid_y = (vertices[:, 1].max() + vertices[:, 1].min()) * 0.5
        mid_z = (vertices[:, 2].max() + vertices[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range * 0.5, mid_x + max_range * 0.5)
        ax.set_ylim(mid_y - max_range * 0.5, mid_y + max_range * 0.5)
        ax.set_zlim(mid_z - max_range * 0.25, mid_z + max_range * 0.25)
        
        # Set title with time information
        time_value = float(glb_path.split('_t')[-1].split('.glb')[0])
        ax.set_title(f'Ocean Surface at Time = {time_value:.1f}s')
        
        # Remove background grid for cleaner look
        ax.grid(False)
        
        # Save the figure
        plt.savefig(output_image_path, bbox_inches='tight')
        plt.close(fig)
        
    except Exception as e:
        print(f"Error rendering {glb_path}: {e}")
        # Create a simple error image
        img = Image.new('RGB', (width, height), color=(255, 0, 0))  # Red for error
        img.save(output_image_path)
    
    return output_image_path

def create_ocean_animation(num_frames=100, time_step=0.1, output_gif="ocean_animation.gif", duration=0.1):
    """Create an animation of ocean waves over time"""
    # Make sure the Rust program is built
    subprocess.run(["cargo", "build"], check=True)
    
    # Create a temporary directory for storing frames
    with tempfile.TemporaryDirectory() as tmp_dir:
        frames = []
        
        # Generate and render each frame
        for i in range(num_frames):
            time_val = i * time_step
            
            # Generate ocean GLB for this time
            glb_path = os.path.join(tmp_dir, f"ocean_t{time_val:.2f}.glb")
            run_ocean_generator(time_val, glb_path)
            
            # Render to image
            img_path = os.path.join(tmp_dir, f"frame_{i:03d}.png")
            render_glb_to_image(glb_path, img_path)
            
            # Read the image for the GIF
            frames.append(imageio.imread(img_path))
            
            print(f"Completed frame {i+1}/{num_frames}")
        
        # Save the animated GIF
        print(f"Creating GIF animation at {output_gif}...")
        imageio.mimsave(output_gif, frames, duration=duration)
        print(f"Animation saved to {output_gif}")

if __name__ == "__main__":
    create_ocean_animation()
    print("Done!")
