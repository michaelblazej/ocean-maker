use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io;
use std::path::PathBuf;

// Import from our library
use ocean_generator::prelude::*;
use ocean_generator::wave;

/// A program to generate tesselated 2D plane-like meshes with trochoidal wave simulation
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Width of the mesh
    #[clap(short, long, default_value = "10.0")]
    width: f32,

    /// length of the mesh
    #[clap(short, long, default_value = "10.0")]
    length: f32,

    /// Number of vertices along width
    #[clap(short = 'x', long, default_value = "10")]
    width_segments: usize,

    /// Number of vertices along length
    #[clap(short = 'y', long, default_value = "10")]
    length_segments: usize,

    /// Apply random perturbation to vertices (0.0 to 1.0)
    #[clap(short, long, default_value = "0.0")]
    noise: f32,

    /// Output file (defaults to stdout if not specified)
    #[clap(short, long)]
    output: Option<PathBuf>,

    /// Output format (obj, stl, raw, or glb)
    #[clap(short, long, default_value = "obj")]
    format: String,
    
    /// Enable trochoidal wave simulation
    #[clap(long)]
    waves: bool,
    
    /// Number of wave components to simulate
    #[clap(long, default_value = "3")]
    wave_count: usize,
    
    /// Wave amplitude (length of waves)
    #[clap(long, default_value = "0.5")]
    amplitude: f32,
    
    /// Wave length (distance between wave crests)
    #[clap(long, default_value = "2.0")]
    wavelength: f32,
    
    /// Steepness of waves (0.0 to 1.0, where 1.0 is maximum steepness)
    #[clap(long, default_value = "0.5")]
    steepness: f32,
    
    /// Wave direction in degrees (0 = along X axis, 90 = along Z axis)
    #[clap(long, default_value = "45.0")]
    direction: f32,
    
    /// Time parameter for the wave animation (in seconds)
    #[clap(long, default_value = "0.0")]
    time: f32,
    
    /// Random seed for wave generation
    #[clap(long, default_value = "42")]
    seed: u64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Create mesh
    let mut mesh = Mesh::new_plane(
        args.width,
        args.length,
        args.width_segments,
        args.length_segments,
        args.noise,
    );
    
    // Apply wave simulation if requested
    if args.waves {
        let wave_params = wave::generate_wave_params(
            args.wave_count,
            args.amplitude,
            args.wavelength,
            args.steepness,
            args.direction,
            args.seed,
        );
        
        println!("Applying trochoidal wave simulation with {} wave components...", args.wave_count);
        mesh.apply_waves(&wave_params, args.time);
    }
    
    // Parse export format
    let format = args.format.parse::<ExportFormat>()
        .map_err(|e| anyhow::anyhow!("Invalid export format: {}", e))?;
    
    // Output to file or stdout
    match args.output {
        Some(path) => {
            let mut file = File::create(&path)
                .with_context(|| format!("Failed to create output file: {}", path.display()))?;
            
            export_mesh(&mesh, format, &mut file)
                .with_context(|| format!("Failed to write to file: {}", path.display()))?;
                
            println!("Mesh exported to: {}", path.display());
        },
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            
            export_mesh(&mesh, format, &mut handle)
                .context("Failed to write to stdout")?;
        }
    }
    
    Ok(())
}
