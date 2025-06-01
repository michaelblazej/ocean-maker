use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

mod tessendorf;

// Use mesh-tools compatibility types
use mesh_tools::compat::{Point3, Vector2, Vector3};

/// Command-line tool to generate an ocean surface using Tessendorf algorithm
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Width of the ocean surface
    #[arg(short, long, default_value_t = 20.0)]
    width: f32,

    /// Length of the ocean surface
    #[arg(short, long, default_value_t = 20.0)]
    length: f32,

    /// Number of segments along the width
    #[arg(long, default_value_t = 128)]
    width_segments: usize,

    /// Number of segments along the length
    #[arg(long, default_value_t = 128)]
    length_segments: usize,

    /// Wind speed for wave generation
    #[arg(long, default_value_t = 5.0)]
    wind_speed: f32,

    /// Wind direction X component
    #[arg(long, default_value_t = 1.0)]
    wind_dir_x: f32,

    /// Wind direction Y component
    #[arg(long, default_value_t = 0.0)]
    wind_dir_y: f32,

    /// Wave amplitude
    #[arg(short, long, default_value_t = 1.0)]
    amplitude: f32,

    /// Wave choppiness (0.0-1.0)
    #[arg(short, long, default_value_t = 0.5)]
    choppiness: f32,

    /// Random seed for wave generation
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Time parameter for wave animation
    #[arg(short, long, default_value_t = 0.0)]
    time: f32,

    /// Output file path
    #[arg(short, long, default_value = "ocean_surface.glb")]
    output: PathBuf,
}

fn main() -> Result<()> {
    // Parse command line arguments
    let args = Args::parse();

    println!("Generating Tessendorf ocean surface...");
    println!("Ocean dimensions: {}x{} with {}x{} segments", 
             args.width, args.length, args.width_segments, args.length_segments);
    println!("Wave parameters: wind speed={}, amplitude={}, choppiness={}, time={}", 
             args.wind_speed, args.amplitude, args.choppiness, args.time);

    // Create ocean parameters
    let params = tessendorf::OceanParams {
        resolution: args.width_segments.max(args.length_segments),
        size: args.width.max(args.length),
        wind_speed: args.wind_speed,
        wind_direction: glam::Vec2::new(args.wind_dir_x, args.wind_dir_y).normalize(),
        amplitude: args.amplitude,
        choppiness: args.choppiness,
        gravity: 9.81,
        depth: 200.0,
        seed: args.seed,
    };

    // Generate displacement field using Tessendorf algorithm
    println!("Computing wave field...");
    let displacement = tessendorf::compute_wave_field(&params, args.time);

    // Create a mesh from the displacement field
    println!("Creating mesh...");
    let (positions, indices, normals, uvs) = create_mesh(
        args.width,
        args.length,
        args.width_segments,
        args.length_segments,
        &displacement,
    );

    // Export the mesh
    println!("Exporting to {}...", args.output.display());
    export_mesh(&positions, &indices, &normals, &uvs, &args.output)?;

    println!("Done!");
    Ok(())
}

/// Create a mesh from the displacement field
fn create_mesh(
    width: f32,
    length: f32,
    width_segments: usize,
    length_segments: usize,
    displacement: &tessendorf::DisplacementField,
) -> (Vec<Point3<f32>>, Vec<u32>, Vec<Vector3<f32>>, Vec<Vector2<f32>>) {
    let grid_x = width_segments;
    let grid_z = length_segments;
    let resolution = displacement.resolution;

    // Calculate grid size
    let segment_width = width / grid_x as f32;
    let segment_length = length / grid_z as f32;

    let half_width = width / 2.0;
    let half_length = length / 2.0;

    // Prepare vertex arrays
    let mut positions = Vec::with_capacity((grid_x + 1) * (grid_z + 1));
    let mut indices = Vec::with_capacity(grid_x * grid_z * 6);
    let mut normals = Vec::with_capacity((grid_x + 1) * (grid_z + 1));
    let mut uvs = Vec::with_capacity((grid_x + 1) * (grid_z + 1));

    // Track min/max height values
    let mut min_height: f32 = f32::MAX;
    let mut max_height: f32 = f32::MIN;

    // Generate vertices
    for iz in 0..=grid_z {
        let z = iz as f32 * segment_length - half_length;
        
        for ix in 0..=grid_x {
            let x = ix as f32 * segment_width - half_width;
            
            // Sample displacement field at this point
            // Map from model coordinates to displacement field coordinates
            let sample_x = ((ix as f32 / grid_x as f32) * (resolution - 1) as f32).round() as usize;
            let sample_z = ((iz as f32 / grid_z as f32) * (resolution - 1) as f32).round() as usize;
            let sample_x = sample_x.min(resolution - 1);
            let sample_z = sample_z.min(resolution - 1);
            
            // Apply displacements from the wave field
            let dx = displacement.displacement_x[sample_x][sample_z];
            let dz = displacement.displacement_y[sample_x][sample_z];
            
            // Get the wave height from displacement field
            let vertex_y = displacement.height[sample_x][sample_z];
            
            // Add position with displacement
            positions.push(mesh_tools::compat::point3::new(x + dx, vertex_y, z + dz));
            
            // Track min/max height
            min_height = min_height.min(vertex_y);
            max_height = max_height.max(vertex_y);
            
            // Add normal
            let nx = displacement.normal_x[sample_x][sample_z];
            let ny = displacement.normal_z[sample_x][sample_z];
            let nz = displacement.normal_y[sample_x][sample_z];
            normals.push(mesh_tools::compat::vector3::new(nx, ny, nz));
            
            // Add texture coordinate
            uvs.push(mesh_tools::compat::vector2::new(
                ix as f32 / grid_x as f32,
                iz as f32 / grid_z as f32,
            ));
        }
    }

    // Generate indices
    for iz in 0..grid_z {
        for ix in 0..grid_x {
            let a = ix + (iz * (grid_x + 1));
            let b = ix + ((iz + 1) * (grid_x + 1));
            let c = (ix + 1) + ((iz + 1) * (grid_x + 1));
            let d = (ix + 1) + (iz * (grid_x + 1));
            
            // First triangle
            indices.push(a as u32);
            indices.push(b as u32);
            indices.push(d as u32);
            
            // Second triangle
            indices.push(b as u32);
            indices.push(c as u32);
            indices.push(d as u32);
        }
    }

    // Log the min/max height values
    println!("Ocean height range: min = {:.4} units, max = {:.4} units, amplitude = {:.4} units", 
             min_height, max_height, (max_height - min_height) / 2.0);

    // Return the mesh data
    (positions, indices, normals, uvs)
}

/// Export the mesh to a GLB file
fn export_mesh(
    positions: &[Point3<f32>],
    indices: &[u32],
    normals: &[Vector3<f32>],
    uvs: &[Vector2<f32>],
    output_path: &PathBuf,
) -> Result<()> {
    use mesh_tools::{GltfBuilder, Triangle};
    
    // Create a new GLTF builder
    let mut builder = GltfBuilder::new();
    
    // Create a blue material for the ocean
    let ocean_material = builder.create_basic_material(
        Some("OceanMaterial".to_string()),
        [0.0, 0.3, 0.8, 1.0], // Blue color
    );
    
    // Convert indices to triangles
    let mut triangles = Vec::with_capacity(indices.len() / 3);
    for i in (0..indices.len()).step_by(3) {
        if i + 2 < indices.len() {
            triangles.push(Triangle::new(
                indices[i], 
                indices[i+1], 
                indices[i+2]
            ));
        }
    }
    
    // Create the mesh using the public API
    let mesh_index = builder.create_simple_mesh(
        Some("OceanSurface".to_string()),
        positions,
        &triangles,
        Some(normals.to_vec()),
        Some(uvs.to_vec()),
        Some(ocean_material)
    );
    
    // Create a node with our mesh
    let node_index = builder.add_node(
        None, 
        Some(mesh_index),
        None,
        None, 
        None
    );
    
    // Add a scene with our node
    builder.add_scene(Some("Ocean Scene".to_string()), Some(vec![node_index]));
    
    // Export to GLB
    builder.export_glb(output_path.to_str().unwrap())?;
    
    println!("Mesh exported to: {}", output_path.display());
    Ok(())
}
