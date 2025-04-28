use glam::{Vec2, Vec3};
use rand::Rng;
use std::collections::HashMap;
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
        let mut rng = rand::thread_rng();
        let mut vertices = Vec::with_capacity((width_segments + 1) * (length_segments + 1));
        let mut faces = Vec::with_capacity(width_segments * length_segments * 2);
        
        // Create vertices
        for y in 0..=length_segments {
            let v = y as f32 / length_segments as f32;
            
            for x in 0..=width_segments {
                let u = x as f32 / width_segments as f32;
                
                // Calculate position
                let mut position = Vec3::new(
                    width * (u - 0.5),
                    0.0,
                    length * (v - 0.5),
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
        for y in 0..length_segments {
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
    pub fn apply_waves(&mut self, wave_params: &[WaveParams], time: f32) {
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
                let phase = proj * (2.0 * std::f32::consts::PI / wave.wavelength) + time * wave.frequency;
                
                // Calculate Gerstner wave displacement
                let steepness = wave.steepness / (2.0 * std::f32::consts::PI / wave.wavelength * wave.amplitude);
                
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
                let nx = -wave.direction.x * wave.amplitude * f32::cos(phase);
                let nz = -wave.direction.y * wave.amplitude * f32::cos(phase);
                let ny = steepness * wave.amplitude * f32::sin(phase);
                
                normal.x -= nx;
                normal.y += ny;
                normal.z -= nz;
            }
            
            // Update the vertex
            vertex.position = final_pos;
            vertex.normal = normal.normalize();
        });
    }
    
    /// Calculate vertex normals based on faces
    pub fn calculate_normals(mut vertices: Vec<Vertex>, faces: &[Face]) -> Vec<Vertex> {
        // Create a map to store accumulated normals
        let mut normal_map: HashMap<usize, Vec3> = HashMap::new();
        
        // Calculate face normals and accumulate them for vertices
        for &Face(i1, i2, i3) in faces {
            let v1 = vertices[i1].position;
            let v2 = vertices[i2].position;
            let v3 = vertices[i3].position;
            
            // Calculate face normal
            let edge1 = v2 - v1;
            let edge2 = v3 - v1;
            let normal = edge1.cross(edge2).normalize();
            
            // Accumulate normals for each vertex
            *normal_map.entry(i1).or_insert(Vec3::ZERO) += normal;
            *normal_map.entry(i2).or_insert(Vec3::ZERO) += normal;
            *normal_map.entry(i3).or_insert(Vec3::ZERO) += normal;
        }
        
        // Apply the calculated normals
        for (i, normal) in normal_map {
            vertices[i].normal = normal.normalize();
        }
        
        vertices
    }
    
    /// Recalculate normals for the existing mesh
    pub fn recalculate_normals(&mut self) {
        // Create a map to store accumulated normals
        let mut normal_map: HashMap<usize, Vec3> = HashMap::new();
        
        // Calculate face normals and accumulate them for vertices
        for &Face(i1, i2, i3) in &self.faces {
            let v1 = self.vertices[i1].position;
            let v2 = self.vertices[i2].position;
            let v3 = self.vertices[i3].position;
            
            // Calculate face normal
            let edge1 = v2 - v1;
            let edge2 = v3 - v1;
            let normal = edge1.cross(edge2).normalize();
            
            // Accumulate normals for each vertex
            *normal_map.entry(i1).or_insert(Vec3::ZERO) += normal;
            *normal_map.entry(i2).or_insert(Vec3::ZERO) += normal;
            *normal_map.entry(i3).or_insert(Vec3::ZERO) += normal;
        }
        
        // Apply the calculated normals
        for (i, normal) in normal_map {
            self.vertices[i].normal = normal.normalize();
        }
    }
}
