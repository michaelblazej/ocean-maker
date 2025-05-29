use std::io::{self, Write};
use std::str::FromStr;
use std::fs::File;
use std::path::Path;
use tempfile::tempdir;

use crate::mesh::Mesh;

/// Only GLB format is supported for export.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExportFormat;

impl FromStr for ExportFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "glb" => Ok(ExportFormat),
            _ => Err(format!("Only GLB format is supported for export, got: {}", s)),
        }
    }
}

/// Export mesh to a file with the specified format.
pub fn export_mesh<W: Write>(mesh: &Mesh, _format: ExportFormat, writer: &mut W) -> io::Result<()> {
    // Only GLB format is supported now
    mesh.export_glb(writer, 1, 1) // Default to 1x1 grid
}

/// Export mesh to a file with the specified format with tiling options.
pub fn export_mesh_tiled<W: Write>(mesh: &Mesh, _format: ExportFormat, writer: &mut W, tile_count_x: usize, tile_count_y: usize) -> io::Result<()> {
    // Only GLB format is supported now
    mesh.export_glb(writer, tile_count_x, tile_count_y)
}

impl Mesh {
    /// Export the mesh to a GLB file
    /// 
    /// # Arguments
    /// * `path` - The path to save the GLB file to
    pub fn save_glb<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut file = File::create(path)?;
        self.export_glb(&mut file, 1, 1) // Default to 1x1 grid
    }
    
    /// Export the mesh to a GLB file with tiling options
    /// 
    /// # Arguments
    /// * `path` - The path to save the GLB file to
    /// * `tile_count_x` - Number of tiles in X direction
    /// * `tile_count_y` - Number of tiles in Y direction
    pub fn save_glb_tiled<P: AsRef<Path>>(&self, path: P, tile_count_x: usize, tile_count_y: usize) -> io::Result<()> {
        let mut file = File::create(path)?;
        self.export_glb(&mut file, tile_count_x, tile_count_y)
    }
    
    /// Export the mesh to a GLB file, writing to the provided writer
    /// 
    /// # Arguments
    /// * `writer` - The writer to write the GLB data to
    /// * `tile_count_x` - Number of tiles in X direction
    /// * `tile_count_y` - Number of tiles in Y direction
    pub fn export_glb<W: Write>(&self, writer: &mut W, tile_count_x: usize, tile_count_y: usize) -> io::Result<()> {
        // Create temporary directory for export
        let temp_dir = tempdir()?;
        let temp_file_path = temp_dir.path().join("ocean.glb");
        
        // Create a new GltfBuilder
        let mut builder = mesh_tools::GltfBuilder::new();
        
        // Create a water material with translucent blue color and metallic reflections
        let water_material = builder.create_metallic_material(
            Some("WaterMaterial".to_string()),
            [0.0, 0.4, 0.8, 0.8],  // Blue water with transparency
            0.9, // Metallic (reflective)
            0.1  // Low roughness (smooth surface)
        );
        
        // Set the material to be double-sided with transparency
        if let Some(materials) = &mut builder.gltf.materials {
            if let Some(material) = materials.get_mut(water_material) {
                material.double_sided = Some(true);
                material.alpha_mode = Some("BLEND".to_string());
            }
        }
        
        // Convert vertex data to format needed by mesh_tools
        let mut positions = Vec::with_capacity(self.vertices.len());
        let mut normals = Vec::with_capacity(self.vertices.len());
        let mut texcoords = Vec::with_capacity(self.vertices.len());
        
        for vertex in &self.vertices {
            // Convert positions
            positions.push(nalgebra::Point3::new(
                vertex.position.x,
                vertex.position.y,
                vertex.position.z
            ));
            
            // Convert normals
            normals.push(nalgebra::Vector3::new(
                vertex.normal.x,
                vertex.normal.y,
                vertex.normal.z
            ));
            
            // Convert UVs
            texcoords.push(nalgebra::Vector2::new(
                vertex.uv.x,
                vertex.uv.y
            ));
        }
        
        // Convert triangle indices
        let mut triangles = Vec::with_capacity(self.faces.len());
        for face in &self.faces {
            triangles.push(mesh_tools::Triangle::new(
                face.0 as u32,
                face.1 as u32,
                face.2 as u32
            ));
        }
        
        // Create a mesh using the create_simple_mesh method
        let mesh_index = builder.create_simple_mesh(
            Some("OceanMesh".to_string()),
            &positions,
            &triangles,
            Some(normals),
            Some(texcoords),
            Some(water_material)
        );
        
        // Create a grid of ocean nodes based on tile counts
        let mut grid_nodes = Vec::new();
        
        // Calculate bounds for placement
        let mut min_x = f32::MAX;
        let mut min_y = f32::MAX;
        let mut max_x = f32::MIN;
        let mut max_y = f32::MIN;
        
        for vertex in &self.vertices {
            min_x = min_x.min(vertex.position.x);
            min_y = min_y.min(vertex.position.y);
            max_x = max_x.max(vertex.position.x);
            max_y = max_y.max(vertex.position.y);
        }
        
        let width = max_x - min_x;
        let depth = max_y - min_y;
        
        // Calculate total grid dimensions
        let grid_width = width * tile_count_x as f32 * 0.9; // Account for 10% overlap
        let grid_height = depth * tile_count_y as f32 * 0.9; // Account for 10% overlap
        
        // Calculate origin offset to center the entire grid
        let origin_x_offset = -grid_width / 2.0 + width * 0.45; // Half grid with half-cell adjustment
        let origin_y_offset = -grid_height / 2.0 + depth * 0.45; // Half grid with half-cell adjustment

        // Create grid of nodes based on provided tile counts
        for row in 0..tile_count_y {
            for col in 0..tile_count_x {
                // Calculate position offset for this grid cell with centering
                let x_offset = origin_x_offset + col as f32 * width * 0.9;  // Slight overlap
                let y_offset = origin_y_offset + row as f32 * depth * 0.9;  // Slight overlap
                
                // Create node for this grid position
                let node_name = format!("OceanNode_{}_{}", row, col);
                let node_index = builder.add_node(
                    Some(node_name),
                    Some(mesh_index),
                    Some([x_offset, y_offset, 0.0]),  // Position at grid coordinates
                    None,  // No rotation
                    None   // No scale
                );
                
                grid_nodes.push(node_index);
            }
        }
        
        // Create a parent node for the grid
        let grid_parent = builder.add_node_with_children(
            Some("OceanGrid".to_string()),
            None,  // No mesh at this level
            None,  // No translation
            None,  // No rotation
            None,  // No scale
            grid_nodes
        );
        
        // Create a scene with the grid parent node
        let scene_index = builder.add_scene(
            Some("OceanScene".to_string()),
            Some(vec![grid_parent])
        );
        
        // Set as default scene
        builder.gltf.scene = Some(scene_index);
        
        // Export to temporary file
        let temp_file_str = temp_file_path.to_str().ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Failed to convert path to string")
        })?;
        
        builder.export_glb(temp_file_str).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("Failed to export GLB: {}", e))
        })?;
        
        // Copy the GLB data to the output writer
        let glb_data = std::fs::read(&temp_file_path)?;
        writer.write_all(&glb_data)?;
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::Mesh;
    use std::fs;
    
    #[test]
    fn test_glb_export() {
        // Create a minimal test mesh (just a single triangle)
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        
        // Create 3 vertices for a single triangle
        vertices.push(crate::mesh::Vertex {
            position: glam::Vec3::new(0.0, 0.0, 0.0),
            normal: glam::Vec3::new(0.0, 0.0, 1.0),
            uv: glam::Vec2::new(0.0, 0.0),
        });
        
        vertices.push(crate::mesh::Vertex {
            position: glam::Vec3::new(1.0, 0.0, 0.0),
            normal: glam::Vec3::new(0.0, 0.0, 1.0),
            uv: glam::Vec2::new(1.0, 0.0),
        });
        
        vertices.push(crate::mesh::Vertex {
            position: glam::Vec3::new(0.0, 1.0, 0.0),
            normal: glam::Vec3::new(0.0, 0.0, 1.0),
            uv: glam::Vec2::new(0.0, 1.0),
        });
        
        // One triangle face
        faces.push(crate::mesh::Face(0, 1, 2));
        
        let mesh = Mesh { vertices, faces };
        
        // Export to a temporary GLB file
        let temp_dir = tempdir().expect("Failed to create temp directory");
        let test_path = temp_dir.path().join("test.glb");
        
        mesh.save_glb(&test_path).expect("Failed to save GLB file");
        
        // Verify the file exists and has content
        assert!(test_path.exists(), "GLB file was not created");
        let metadata = fs::metadata(&test_path).expect("Failed to get file metadata");
        assert!(metadata.len() > 0, "GLB file is empty");
    }
}
