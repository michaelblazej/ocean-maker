use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;

use ocean_generator::prelude::*;
use ocean_generator::wave;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Generating animated ocean surface...");
    
    // Parameters for the ocean surface
    let width = 20.0;
    let length = 20.0;
    let resolution = 50; // Higher resolution for smoother surface
    
    // Wave parameters
    let wave_count = 5;
    let amplitude = 0.7;
    let wavelength = 4.0;
    let steepness = 0.6;
    let direction = 45.0; // 45 degrees
    let seed = 42;
    
    // Create base mesh
    let mut mesh = Mesh::new_plane(
        width,
        length,
        resolution,
        resolution,
        0.0, // No noise
    );
    
    // Generate wave parameters
    let wave_params = wave::generate_wave_params(
        wave_count,
        amplitude,
        wavelength,
        steepness,
        direction,
        seed,
    );
    
    // Animation parameters
    let frame_count = 60;  // Generate 60 frames (one second at 60fps)
    let time_step = 0.05;  // Time step between frames
    
    // Create output directory
    let output_dir = PathBuf::from("./animation_output");
    std::fs::create_dir_all(&output_dir)?;
    
    println!("Generating {} frames of animation...", frame_count);
    let start_time = Instant::now();
    
    // Generate frames
    for frame in 0..frame_count {
        // Calculate animation time
        let time = frame as f32 * time_step;
        
        // Create a copy of the base mesh for this frame
        let mut frame_mesh = mesh.clone();
        
        // Apply wave simulation at this time
        frame_mesh.apply_waves(&wave_params, time);
        
        // Export to OBJ file
        let filename = format!("ocean_frame_{:03}.obj", frame);
        let file_path = output_dir.join(filename);
        
        let mut file = File::create(&file_path)?;
        export_mesh(&frame_mesh, ExportFormat::Obj, &mut file)?;
        
        if frame % 10 == 0 {
            println!("Generated frame {}/{}", frame + 1, frame_count);
        }
    }
    
    let elapsed = start_time.elapsed();
    println!("Animation generation complete in {:.2?}", elapsed);
    println!("Output files saved to: {}", output_dir.display());
    
    Ok(())
}
