use std::f32::consts::PI;
use std::fs::File;
use std::path::PathBuf;

use ocean_generator::prelude::*;
use ocean_generator::wave;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Generating ocean with floating object...");
    
    // Create ocean mesh
    let mut ocean = Mesh::new_plane(
        30.0,       // width
        30.0,       // length
        60,         // width segments
        60,         // length segments
        0.0,        // no random noise
    );
    
    // Create simple boat mesh (just a rectangular prism)
    let boat = create_boat_mesh(3.0, 1.5, 0.8);
    
    // Generate wave parameters for the ocean
    let wave_params = wave::generate_wave_params(
        4,          // 4 wave components
        0.5,        // amplitude
        6.0,        // wavelength
        0.4,        // steepness
        30.0,       // direction (degrees)
        123,        // random seed
    );
    
    // Apply wave simulation to the ocean at time=0
    ocean.apply_waves(&wave_params, 0.0);
    
    // Sample points on the ocean surface to determine boat position and orientation
    let time = 0.0;
    let boat_center_x = 5.0;
    let boat_center_z = 2.0;
    
    // Create output directory
    let output_dir = PathBuf::from("./floating_output");
    std::fs::create_dir_all(&output_dir)?;
    
    // Export ocean
    let ocean_path = output_dir.join("ocean.obj");
    let mut ocean_file = File::create(&ocean_path)?;
    export_mesh(&ocean, ExportFormat::Obj, &mut ocean_file)?;
    
    // Export boat positioned on the ocean
    let boat_path = output_dir.join("boat.obj");
    let mut boat_file = File::create(&boat_path)?;
    export_mesh(&boat, ExportFormat::Obj, &mut boat_file)?;
    
    println!("Ocean and floating object exported to: {}", output_dir.display());
    
    Ok(())
}

fn create_boat_mesh(length: f32, width: f32, height: f32) -> Mesh {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    
    // Create a simple boat shape (rectangular prism)
    let half_length = length / 2.0;
    let half_width = width / 2.0;
    
    // Bottom vertices
    vertices.push(Vertex {
        position: glam::Vec3::new(-half_length, 0.0, -half_width),
        normal: glam::Vec3::new(0.0, -1.0, 0.0),
        uv: glam::Vec2::new(0.0, 0.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(half_length, 0.0, -half_width),
        normal: glam::Vec3::new(0.0, -1.0, 0.0),
        uv: glam::Vec2::new(1.0, 0.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(half_length, 0.0, half_width),
        normal: glam::Vec3::new(0.0, -1.0, 0.0),
        uv: glam::Vec2::new(1.0, 1.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(-half_length, 0.0, half_width),
        normal: glam::Vec3::new(0.0, -1.0, 0.0),
        uv: glam::Vec2::new(0.0, 1.0),
    });
    
    // Top vertices
    vertices.push(Vertex {
        position: glam::Vec3::new(-half_length, height, -half_width),
        normal: glam::Vec3::new(0.0, 1.0, 0.0),
        uv: glam::Vec2::new(0.0, 0.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(half_length, height, -half_width),
        normal: glam::Vec3::new(0.0, 1.0, 0.0),
        uv: glam::Vec2::new(1.0, 0.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(half_length, height, half_width),
        normal: glam::Vec3::new(0.0, 1.0, 0.0),
        uv: glam::Vec2::new(1.0, 1.0),
    });
    vertices.push(Vertex {
        position: glam::Vec3::new(-half_length, height, half_width),
        normal: glam::Vec3::new(0.0, 1.0, 0.0),
        uv: glam::Vec2::new(0.0, 1.0),
    });
    
    // Bottom face
    faces.push(Face(0, 2, 1));
    faces.push(Face(0, 3, 2));
    
    // Top face
    faces.push(Face(4, 5, 6));
    faces.push(Face(4, 6, 7));
    
    // Side faces
    faces.push(Face(0, 1, 5));
    faces.push(Face(0, 5, 4));
    
    faces.push(Face(1, 2, 6));
    faces.push(Face(1, 6, 5));
    
    faces.push(Face(2, 3, 7));
    faces.push(Face(2, 7, 6));
    
    faces.push(Face(3, 0, 4));
    faces.push(Face(3, 4, 7));
    
    Mesh { vertices, faces }
}
