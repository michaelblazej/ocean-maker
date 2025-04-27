use anyhow::{Context, Result};
use byteorder::{LittleEndian, WriteBytesExt};
use clap::Parser;
use glam::{Vec2, Vec3};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rayon::prelude::*;
use serde_json::json;
use std::f32::consts::PI;
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

/// A program to generate tesselated 2D plane-like meshes with trochoidal wave simulation
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Width of the mesh
    #[clap(short, long, default_value = "10.0")]
    width: f32,

    /// Height of the mesh
    #[clap(short = 'H', long, default_value = "10.0")]
    height: f32,

    /// Number of vertices along width
    #[clap(short = 'x', long, default_value = "10")]
    width_segments: usize,

    /// Number of vertices along height
    #[clap(short = 'y', long, default_value = "10")]
    height_segments: usize,

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
    
    /// Wave amplitude (height of waves)
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

/// A vertex in 3D space
#[derive(Debug, Clone, Copy)]
struct Vertex {
    position: Vec3,
    normal: Vec3,
    uv: Vec2,
}

/// A face consisting of three vertex indices
#[derive(Debug, Clone, Copy)]
struct Face(usize, usize, usize);

/// Parameters for a single wave component
#[derive(Debug, Clone, Copy)]
struct WaveParams {
    amplitude: f32,
    wavelength: f32,
    steepness: f32,
    direction: Vec2,
    frequency: f32,
}

/// The complete mesh
struct Mesh {
    vertices: Vec<Vertex>,
    faces: Vec<Face>,
}

impl Mesh {
    /// Create a new plane mesh with the given dimensions and resolution
    fn new_plane(width: f32, height: f32, width_segments: usize, height_segments: usize, noise: f32) -> Self {
        let mut rng = rand::thread_rng();
        let mut vertices = Vec::with_capacity((width_segments + 1) * (height_segments + 1));
        let mut faces = Vec::with_capacity(width_segments * height_segments * 2);
        
        // Create vertices
        for y in 0..=height_segments {
            let v = y as f32 / height_segments as f32;
            
            for x in 0..=width_segments {
                let u = x as f32 / width_segments as f32;
                
                // Calculate position
                let mut position = Vec3::new(
                    width * (u - 0.5),
                    0.0,
                    height * (v - 0.5),
                );
                
                // Apply noise if requested
                if noise > 0.0 {
                    position.y = (rng.gen::<f32>() - 0.5) * noise;
                }
                
                vertices.push(Vertex {
                    position,
                    normal: Vec3::new(0.0, 1.0, 0.0),
                    uv: Vec2::new(u, v),
                });
            }
        }
        
        // Create faces (triangles)
        for y in 0..height_segments {
            for x in 0..width_segments {
                let a = x + y * (width_segments + 1);
                let b = x + (y + 1) * (width_segments + 1);
                let c = (x + 1) + (y + 1) * (width_segments + 1);
                let d = (x + 1) + y * (width_segments + 1);
                
                // First triangle
                faces.push(Face(a, b, d));
                // Second triangle
                faces.push(Face(b, c, d));
            }
        }
        
        // Recalculate normals if we applied noise
        let vertices = if noise > 0.0 {
            Self::calculate_normals(vertices, &faces)
        } else {
            vertices
        };
        
        Mesh { vertices, faces }
    }
    
    /// Apply trochoidal wave simulation to the mesh
    fn apply_waves(&mut self, wave_params: &[WaveParams], time: f32) {
        // Apply wave deformation to each vertex in parallel
        self.vertices.par_iter_mut().for_each(|vertex| {
            let pos = vertex.position;
            // Start with the original position for accumulation
            let mut final_pos = pos;
            let mut normal = Vec3::new(0.0, 1.0, 0.0);
            
            for wave in wave_params {
                // Project position onto wave direction
                let pos_2d = Vec2::new(pos.x, pos.z);
                let proj = wave.direction.dot(pos_2d);
                
                // Calculate phase
                let phase = proj * (2.0 * PI / wave.wavelength) + time * wave.frequency;
                
                // Calculate Gerstner wave displacement
                let steepness = wave.steepness / (2.0 * PI / wave.wavelength * wave.amplitude);
                
                // Horizontal displacement
                let dx = steepness * wave.amplitude * wave.direction.x * f32::cos(phase);
                let dz = steepness * wave.amplitude * wave.direction.y * f32::cos(phase);
                
                // Vertical displacement
                let dy = wave.amplitude * f32::sin(phase);
                
                // Accumulate displacement
                final_pos.x += dx;
                final_pos.y += dy;
                final_pos.z += dz;
                
                // Accumulate normal contribution
                let normal_x = -wave.direction.x * wave.steepness * f32::cos(phase);
                let normal_z = -wave.direction.y * wave.steepness * f32::cos(phase);
                normal.x -= normal_x;
                normal.z -= normal_z;
            }
            
            // Update vertex position
            vertex.position = final_pos;
            
            // Normalize and update normal
            normal = normal.normalize();
            vertex.normal = normal;
        });
        
        // Calculate new normals based on face geometry for better results
        self.recalculate_normals();
    }
    
    /// Generate random wave parameters
    fn generate_wave_params(
        count: usize,
        base_amplitude: f32,
        base_wavelength: f32,
        base_steepness: f32,
        base_direction_degrees: f32,
        seed: u64,
    ) -> Vec<WaveParams> {
        let mut rng: StdRng = SeedableRng::seed_from_u64(seed);
        let base_direction_radians = base_direction_degrees * PI / 180.0;
        
        let mut waves = Vec::with_capacity(count);
        
        for i in 0..count {
            // Each successive wave has:
            // - Smaller amplitude (1/i+1 falloff)
            // - Different wavelength (with some randomness)
            // - Varied direction (within ±30 degrees of base direction)
            // - Randomized phase
            
            let amplitude_factor = 1.0 / (i as f32 + 1.0);
            let amplitude = base_amplitude * amplitude_factor;
            
            // Wavelength varies between 0.5x and 1.5x base wavelength
            let wavelength = base_wavelength * (0.5 + rng.gen::<f32>());
            
            // Steepness is consistent but slightly randomized
            let steepness = base_steepness * (0.8 + 0.4 * rng.gen::<f32>());
            
            // Direction varies within ±30 degrees of base
            let angle_offset = (rng.gen::<f32>() - 0.5) * 60.0 * PI / 180.0;
            let direction_angle = base_direction_radians + angle_offset;
            let direction = Vec2::new(direction_angle.cos(), direction_angle.sin());
            
            // Calculate frequency (for time animation)
            let frequency = (2.0 * PI / wavelength).sqrt() * 1.5;
            
            waves.push(WaveParams {
                amplitude,
                wavelength,
                steepness,
                direction,
                frequency,
            });
        }
        
        waves
    }
    
    /// Calculate vertex normals based on faces
    fn calculate_normals(mut vertices: Vec<Vertex>, faces: &[Face]) -> Vec<Vertex> {
        // Reset normals
        for vertex in &mut vertices {
            vertex.normal = Vec3::ZERO;
        }
        
        // Calculate face normals and add to vertices
        for face in faces {
            let v0 = vertices[face.0].position;
            let v1 = vertices[face.1].position;
            let v2 = vertices[face.2].position;
            
            let normal = (v1 - v0).cross(v2 - v0).normalize();
            
            vertices[face.0].normal += normal;
            vertices[face.1].normal += normal;
            vertices[face.2].normal += normal;
        }
        
        // Normalize all vertex normals
        for vertex in &mut vertices {
            if vertex.normal != Vec3::ZERO {
                vertex.normal = vertex.normal.normalize();
            } else {
                vertex.normal = Vec3::Y;
            }
        }
        
        vertices
    }
    
    /// Recalculate normals for the existing mesh
    fn recalculate_normals(&mut self) {
        // Reset normals
        for vertex in &mut self.vertices {
            vertex.normal = Vec3::ZERO;
        }
        
        // Calculate face normals and add to vertices
        for face in &self.faces {
            let v0 = self.vertices[face.0].position;
            let v1 = self.vertices[face.1].position;
            let v2 = self.vertices[face.2].position;
            
            let normal = (v1 - v0).cross(v2 - v0).normalize();
            
            self.vertices[face.0].normal += normal;
            self.vertices[face.1].normal += normal;
            self.vertices[face.2].normal += normal;
        }
        
        // Normalize all vertex normals
        for vertex in &mut self.vertices {
            if vertex.normal != Vec3::ZERO {
                vertex.normal = vertex.normal.normalize();
            } else {
                vertex.normal = Vec3::Y;
            }
        }
    }
    
    /// Export mesh to OBJ format
    fn export_obj<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Write header
        writeln!(writer, "# Generated by ocean-generator")?;
        
        // Write vertices
        for vertex in &self.vertices {
            writeln!(
                writer,
                "v {:.6} {:.6} {:.6}",
                vertex.position.x, vertex.position.y, vertex.position.z
            )?;
        }
        
        // Write texture coordinates
        for vertex in &self.vertices {
            writeln!(writer, "vt {:.6} {:.6}", vertex.uv.x, vertex.uv.y)?;
        }
        
        // Write normals
        for vertex in &self.vertices {
            writeln!(
                writer,
                "vn {:.6} {:.6} {:.6}",
                vertex.normal.x, vertex.normal.y, vertex.normal.z
            )?;
        }
        
        // Write faces
        for face in &self.faces {
            writeln!(
                writer,
                "f {}/{}/{} {}/{}/{} {}/{}/{}",
                face.0 + 1, face.0 + 1, face.0 + 1,
                face.1 + 1, face.1 + 1, face.1 + 1,
                face.2 + 1, face.2 + 1, face.2 + 1
            )?;
        }
        
        Ok(())
    }
    
    /// Export mesh to STL format
    fn export_stl<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Write header
        writeln!(writer, "solid generated_mesh")?;
        
        for face in &self.faces {
            let v0 = self.vertices[face.0].position;
            let v1 = self.vertices[face.1].position;
            let v2 = self.vertices[face.2].position;
            
            let normal = self.vertices[face.0].normal;
            
            writeln!(writer, "  facet normal {:.6} {:.6} {:.6}", normal.x, normal.y, normal.z)?;
            writeln!(writer, "    outer loop")?;
            writeln!(writer, "      vertex {:.6} {:.6} {:.6}", v0.x, v0.y, v0.z)?;
            writeln!(writer, "      vertex {:.6} {:.6} {:.6}", v1.x, v1.y, v1.z)?;
            writeln!(writer, "      vertex {:.6} {:.6} {:.6}", v2.x, v2.y, v2.z)?;
            writeln!(writer, "    endloop")?;
            writeln!(writer, "  endfacet")?;
        }
        
        writeln!(writer, "endsolid generated_mesh")?;
        Ok(())
    }
    
    /// Export mesh to a simple raw format (vertices and indices)
    fn export_raw<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Write header with counts
        writeln!(writer, "vertices: {}", self.vertices.len())?;
        writeln!(writer, "faces: {}", self.faces.len())?;
        
        // Write vertices
        for vertex in &self.vertices {
            writeln!(
                writer,
                "v {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
                vertex.position.x, vertex.position.y, vertex.position.z,
                vertex.normal.x, vertex.normal.y, vertex.normal.z,
                vertex.uv.x, vertex.uv.y
            )?;
        }
        
        // Write faces
        for face in &self.faces {
            writeln!(writer, "f {} {} {}", face.0, face.1, face.2)?;
        }
        
        Ok(())
    }
    
    /// Export mesh to GLB format (GL Transmission Format Binary)
    fn export_glb<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Create buffers for vertex and index data
        let mut vertex_buffer = Vec::new();
        let mut normal_buffer = Vec::new();
        let mut texcoord_buffer = Vec::new();
        let mut index_buffer = Vec::new();
        
        // Populate vertex buffers
        for vertex in &self.vertices {
            // Position (convert to right-handed coordinate system for glTF)
            vertex_buffer.push(vertex.position.x);
            vertex_buffer.push(vertex.position.y);
            vertex_buffer.push(vertex.position.z);
            
            // Normal
            normal_buffer.push(vertex.normal.x);
            normal_buffer.push(vertex.normal.y);
            normal_buffer.push(vertex.normal.z);
            
            // Texture coordinates
            texcoord_buffer.push(vertex.uv.x);
            texcoord_buffer.push(vertex.uv.y);
        }
        
        // Populate index buffer
        for face in &self.faces {
            index_buffer.push(face.0 as u32);
            index_buffer.push(face.1 as u32);
            index_buffer.push(face.2 as u32);
        }
        
        // Create binary buffer for vertex and index data
        let mut buffer_data = Vec::new();
        
        // Add vertex positions
        for value in &vertex_buffer {
            buffer_data.write_f32::<LittleEndian>(*value)?;
        }
        
        // Add normals
        for value in &normal_buffer {
            buffer_data.write_f32::<LittleEndian>(*value)?;
        }
        
        // Add texture coordinates
        for value in &texcoord_buffer {
            buffer_data.write_f32::<LittleEndian>(*value)?;
        }
        
        // Add indices
        for value in &index_buffer {
            buffer_data.write_u32::<LittleEndian>(*value)?;
        }
        
        // Pad the buffer to 4-byte alignment
        while buffer_data.len() % 4 != 0 {
            buffer_data.push(0);
        }
        
        // Create JSON for glTF
        let vertex_count = self.vertices.len();
        let face_count = self.faces.len();
        let positions_byte_offset = 0;
        let normals_byte_offset = vertex_count * 3 * std::mem::size_of::<f32>();
        let texcoords_byte_offset = normals_byte_offset + vertex_count * 3 * std::mem::size_of::<f32>();
        let indices_byte_offset = texcoords_byte_offset + vertex_count * 2 * std::mem::size_of::<f32>();
        
        let json = json!({
            "asset": {
                "version": "2.0",
                "generator": "ocean-generator"
            },
            "scene": 0,
            "scenes": [{
                "nodes": [0]
            }],
            "nodes": [{
                "mesh": 0
            }],
            "meshes": [{
                "primitives": [{
                    "attributes": {
                        "POSITION": 0,
                        "NORMAL": 1,
                        "TEXCOORD_0": 2
                    },
                    "indices": 3,
                    "mode": 4  // TRIANGLES
                }]
            }],
            "buffers": [{
                "byteLength": buffer_data.len()
            }],
            "bufferViews": [
                {
                    "buffer": 0,
                    "byteOffset": positions_byte_offset,
                    "byteLength": vertex_count * 3 * std::mem::size_of::<f32>(),
                    "target": 34962  // ARRAY_BUFFER
                },
                {
                    "buffer": 0,
                    "byteOffset": normals_byte_offset,
                    "byteLength": vertex_count * 3 * std::mem::size_of::<f32>(),
                    "target": 34962  // ARRAY_BUFFER
                },
                {
                    "buffer": 0,
                    "byteOffset": texcoords_byte_offset,
                    "byteLength": vertex_count * 2 * std::mem::size_of::<f32>(),
                    "target": 34962  // ARRAY_BUFFER
                },
                {
                    "buffer": 0,
                    "byteOffset": indices_byte_offset,
                    "byteLength": face_count * 3 * std::mem::size_of::<u32>(),
                    "target": 34963  // ELEMENT_ARRAY_BUFFER
                }
            ],
            "accessors": [
                {
                    "bufferView": 0,
                    "componentType": 5126,  // FLOAT
                    "count": vertex_count,
                    "type": "VEC3",
                    "max": [
                        vertex_buffer.chunks(3).map(|v| v[0]).fold(f32::NEG_INFINITY, f32::max),
                        vertex_buffer.chunks(3).map(|v| v[1]).fold(f32::NEG_INFINITY, f32::max),
                        vertex_buffer.chunks(3).map(|v| v[2]).fold(f32::NEG_INFINITY, f32::max)
                    ],
                    "min": [
                        vertex_buffer.chunks(3).map(|v| v[0]).fold(f32::INFINITY, f32::min),
                        vertex_buffer.chunks(3).map(|v| v[1]).fold(f32::INFINITY, f32::min),
                        vertex_buffer.chunks(3).map(|v| v[2]).fold(f32::INFINITY, f32::min)
                    ]
                },
                {
                    "bufferView": 1,
                    "componentType": 5126,  // FLOAT
                    "count": vertex_count,
                    "type": "VEC3"
                },
                {
                    "bufferView": 2,
                    "componentType": 5126,  // FLOAT
                    "count": vertex_count,
                    "type": "VEC2"
                },
                {
                    "bufferView": 3,
                    "componentType": 5125,  // UNSIGNED_INT
                    "count": face_count * 3,
                    "type": "SCALAR"
                }
            ]
        });
        
        // Convert JSON to string and pad to 4-byte alignment
        let mut json_string = serde_json::to_string(&json)?;
        while json_string.len() % 4 != 0 {
            json_string.push(' ');
        }
        
        // Write GLB header
        let glb_length = 12 + 8 + json_string.len() + 8 + buffer_data.len();
        writer.write_u32::<LittleEndian>(0x46546C67)?; // magic: glTF
        writer.write_u32::<LittleEndian>(2)?; // version
        writer.write_u32::<LittleEndian>(glb_length as u32)?; // length
        
        // Write JSON chunk
        writer.write_u32::<LittleEndian>(json_string.len() as u32)?; // chunk length
        writer.write_u32::<LittleEndian>(0x4E4F534A)?; // chunk type: JSON
        writer.write_all(json_string.as_bytes())?;
        
        // Write binary chunk
        writer.write_u32::<LittleEndian>(buffer_data.len() as u32)?; // chunk length
        writer.write_u32::<LittleEndian>(0x004E4942)?; // chunk type: BIN
        writer.write_all(&buffer_data)?;
        
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    // Validate parameters
    if args.width <= 0.0 || args.height <= 0.0 {
        anyhow::bail!("Width and height must be positive");
    }
    
    if args.width_segments == 0 || args.height_segments == 0 {
        anyhow::bail!("Width and height segments must be at least 1");
    }
    
    if args.noise < 0.0 || args.noise > 1.0 {
        anyhow::bail!("Noise must be between 0.0 and 1.0");
    }
    
    if args.waves {
        if args.steepness < 0.0 || args.steepness > 1.0 {
            anyhow::bail!("Wave steepness must be between 0.0 and 1.0");
        }
        
        if args.wave_count == 0 {
            anyhow::bail!("Wave count must be at least 1");
        }
        
        if args.amplitude <= 0.0 {
            anyhow::bail!("Wave amplitude must be positive");
        }
        
        if args.wavelength <= 0.0 {
            anyhow::bail!("Wave length must be positive");
        }
    }
    
    // Generate the base mesh
    let mut mesh = Mesh::new_plane(
        args.width,
        args.height,
        args.width_segments,
        args.height_segments,
        args.noise,
    );
    
    // Apply trochoidal wave simulation if enabled
    if args.waves {
        let wave_params = Mesh::generate_wave_params(
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
    
    // Create output writer
    match args.output {
        Some(path) => {
            let mut file = File::create(&path)
                .with_context(|| format!("Failed to create output file: {}", path.display()))?;
            
            export_mesh(&mesh, &args.format, &mut file)
                .with_context(|| format!("Failed to write to file: {}", path.display()))?;
                
            println!("Mesh exported to: {}", path.display());
        },
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            
            export_mesh(&mesh, &args.format, &mut handle)
                .context("Failed to write to stdout")?;
        }
    }
    
    Ok(())
}

fn export_mesh<W: Write>(mesh: &Mesh, format: &str, writer: &mut W) -> io::Result<()> {
    match format.to_lowercase().as_str() {
        "obj" => mesh.export_obj(writer),
        "stl" => mesh.export_stl(writer),
        "raw" => mesh.export_raw(writer),
        "glb" => mesh.export_glb(writer),
        _ => {
            writeln!(writer, "Error: Unsupported format '{}'", format)?;
            writeln!(writer, "Supported formats: obj, stl, raw, glb")?;
            Err(io::Error::new(io::ErrorKind::InvalidInput, "Unsupported format"))
        }
    }
}
