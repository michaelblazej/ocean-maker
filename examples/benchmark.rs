use std::time::Instant;
use ocean_generator::prelude::*;
use ocean_generator::wave;

fn main() {
    println!("Running Ocean Generator Benchmarks");
    println!("==================================");
    
    // Parameters for different mesh sizes
    let resolutions = [
        (100, 100),   // 10k vertices
        (500, 500),   // 250k vertices
        (1000, 1000), // 1M vertices
    ];
    
    for &(width_segments, length_segments) in &resolutions {
        println!("\nMesh size: {}x{} ({} vertices)", 
            width_segments, 
            length_segments,
            (width_segments + 1) * (length_segments + 1)
        );
        
        // Benchmark mesh generation
        let start = Instant::now();
        let mesh = Mesh::new_plane(
            100.0, 
            100.0, 
            width_segments, 
            length_segments, 
            0.1 // slight noise to trigger normal calculation
        );
        let elapsed = start.elapsed();
        println!("  Mesh generation: {:.2?}", elapsed);
        
        // Benchmark wave simulation
        let wave_params = wave::generate_wave_params(
            5,       // 5 wave components
            0.5,     // amplitude
            10.0,    // wavelength
            0.5,     // steepness
            45.0,    // direction
            42,      // seed
        );
        
        let mut wave_mesh = mesh.clone();
        let start = Instant::now();
        wave_mesh.apply_waves(&wave_params, 0.0);
        let elapsed = start.elapsed();
        println!("  Wave simulation: {:.2?}", elapsed);
        
        // Only run export benchmark on smaller meshes to avoid excessive memory usage
        if width_segments <= 500 {
            // Benchmark export
            let start = Instant::now();
            let mut buffer = Vec::new();
            export_mesh(&wave_mesh, ExportFormat::Obj, &mut buffer).unwrap();
            let elapsed = start.elapsed();
            println!("  OBJ export: {:.2?} ({}MB)", elapsed, buffer.len() / (1024 * 1024));
        }
    }
}
