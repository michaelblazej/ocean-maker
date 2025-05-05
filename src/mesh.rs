use glam::{Vec2, Vec3};
use rand::Rng;
use rayon::prelude::*;

use crate::wave::WaveParams;

/// A vertex in 3D space
#[derive(Debug, Clone, Copy)]
pub struct Vertex {
    pub position: Vec3,
    pub normal: Vec3,
    pub uv: Vec2,
}

/// A face consisting of three vertex indices
#[derive(Debug, Clone, Copy)]
pub struct Face(pub usize, pub usize, pub usize);

/// The complete mesh
#[derive(Debug, Clone)]
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub faces: Vec<Face>,
}

impl Mesh {
    /// Create a new plane mesh with the given dimensions and resolution
    pub fn new_plane(width: f32, length: f32, width_segments: usize, length_segments: usize, noise: f32) -> Self {
        let _vertex_count = (width_segments + 1) * (length_segments + 1);
        let _face_count = width_segments * length_segments * 2;
        
        // Create vertices in parallel
        let vertices: Vec<Vertex> = (0..=length_segments)
            .into_par_iter()
            .flat_map(|y| {
                let v = y as f32 / length_segments as f32;
                let mut thread_rng = rand::thread_rng();
                
                (0..=width_segments).map(move |x| {
                    let u = x as f32 / width_segments as f32;
                    
                    // Calculate position
                    let mut position = Vec3::new(
                        width * (u - 0.5),
                        length * (v - 0.5),
                        0.0,
                    );
                    
                    // Apply noise if requested
                    if noise > 0.0 {
                        position.z = (thread_rng.gen::<f32>() - 0.5) * noise;
                    }
                    
                    Vertex {
                        position,
                        normal: Vec3::new(0.0, 0.0, 1.0),
                        uv: Vec2::new(u, v),
                    }
                }).collect::<Vec<_>>()
            })
            .collect();
        
        // Create face indices
        let faces: Vec<Face> = (0..length_segments)
            .into_par_iter()
            .flat_map(|y| {
                (0..width_segments).flat_map(move |x| {
                    let a = x + y * (width_segments + 1);
                    let b = x + (y + 1) * (width_segments + 1);
                    let c = (x + 1) + (y + 1) * (width_segments + 1);
                    let d = (x + 1) + y * (width_segments + 1);
                    
                    // Return both triangles
                    vec![
                        Face(a, b, d),
                        Face(b, c, d)
                    ]
                }).collect::<Vec<_>>()
            })
            .collect();
        
        // Recalculate normals if we applied noise
        let vertices = if noise > 0.0 {
            Self::calculate_normals_parallel(vertices, &faces)
        } else {
            vertices
        };
        
        Mesh { vertices, faces }
    }
    
    /// Calculate vertex normals based on faces (parallel version)
    pub fn calculate_normals_parallel(vertices: Vec<Vertex>, faces: &[Face]) -> Vec<Vertex> {
        let face_normals: Vec<(usize, usize, usize, Vec3)> = faces.par_iter().map(|&Face(i1, i2, i3)| {
            let v1 = vertices[i1].position;
            let v2 = vertices[i2].position;
            let v3 = vertices[i3].position;
            
            // Calculate face normal
            let edge1 = v2 - v1;
            let edge2 = v3 - v1;
            let normal = edge1.cross(edge2).normalize();
            
            (i1, i2, i3, normal)
        }).collect();
        
        // Create a thread-safe accumulator for normals
        let normal_accumulators = std::sync::Mutex::new(vec![Vec3::ZERO; vertices.len()]);
        
        // Accumulate normals in parallel
        face_normals.par_iter().for_each(|&(i1, i2, i3, normal)| {
            let mut accumulators = normal_accumulators.lock().unwrap();
            accumulators[i1] += normal;
            accumulators[i2] += normal;
            accumulators[i3] += normal;
        });
        
        // Apply the accumulated normals to create the final vertices
        let normal_accumulators = normal_accumulators.into_inner().unwrap();
        
        vertices.into_iter().enumerate().map(|(i, mut vertex)| {
            vertex.normal = normal_accumulators[i].normalize();
            vertex
        }).collect()
    }
    
    /// Recalculate normals for the existing mesh (parallel version)
    pub fn recalculate_normals(&mut self) {
        let face_normals: Vec<(usize, usize, usize, Vec3)> = self.faces.par_iter().map(|&Face(i1, i2, i3)| {
            let v1 = self.vertices[i1].position;
            let v2 = self.vertices[i2].position;
            let v3 = self.vertices[i3].position;
            
            // Calculate face normal
            let edge1 = v2 - v1;
            let edge2 = v3 - v1;
            let normal = edge1.cross(edge2).normalize();
            
            (i1, i2, i3, normal)
        }).collect();
        
        // Create a thread-safe accumulator for normals
        let normal_accumulators = std::sync::Mutex::new(vec![Vec3::ZERO; self.vertices.len()]);
        
        // Accumulate normals in parallel
        face_normals.par_iter().for_each(|&(i1, i2, i3, normal)| {
            let mut accumulators = normal_accumulators.lock().unwrap();
            accumulators[i1] += normal;
            accumulators[i2] += normal;
            accumulators[i3] += normal;
        });
        
        // Apply the accumulated normals
        let normal_accumulators = normal_accumulators.into_inner().unwrap();
        
        self.vertices.par_iter_mut().enumerate().for_each(|(i, vertex)| {
            vertex.normal = normal_accumulators[i].normalize();
        });
    }
    
    /// Apply trochoidal wave simulation to the mesh
    pub fn apply_waves(&mut self, wave_params: &[WaveParams], time: f32) {
        // Apply wave deformation to each vertex in parallel
        self.vertices.par_iter_mut().for_each(|vertex| {
            let pos = vertex.position;
            // Start with the original position for accumulation
            let mut final_pos = pos;
            let mut normal = Vec3::new(0.0, 0.0, 1.0);
            
            for wave in wave_params {
                // Project position onto wave direction
                let pos_2d = Vec2::new(pos.x, pos.y);
                let proj = wave.direction.dot(pos_2d);
                
                // Calculate phase
                let phase = proj * (2.0 * std::f32::consts::PI / wave.wavelength) + time * wave.frequency;
                
                // Calculate Gerstner wave displacement
                let steepness = wave.steepness / (2.0 * std::f32::consts::PI / wave.wavelength * wave.amplitude);
                
                // Horizontal displacement
                let dx = steepness * wave.amplitude * wave.direction.x * f32::cos(phase);
                let dy = steepness * wave.amplitude * wave.direction.y * f32::sin(phase);
                
                // Vertical displacement
                let dz = wave.amplitude * f32::sin(phase);
                
                // Accumulate displacement
                final_pos.x += dx;
                final_pos.y += dy;
                final_pos.z += dz;
                
                // Accumulate normal contribution
                let nx = -wave.direction.x * wave.amplitude * f32::cos(phase);
                let ny = -wave.direction.y * wave.amplitude * f32::cos(phase);
                let nz = steepness * wave.amplitude * f32::sin(phase);
                
                normal.x -= nx;
                normal.y += ny;
                normal.z -= nz;
            }
            
            // Update the vertex
            vertex.position = final_pos;
            vertex.normal = normal.normalize();
        });
    }
}
