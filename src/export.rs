use std::io::{self, Write};
use byteorder::{LittleEndian, WriteBytesExt};
use gltf::json;
use std::str::FromStr;

use crate::mesh::Mesh;

/// Supported export formats
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExportFormat {
    Obj,
    Stl,
    Raw,
    Glb,
}

impl FromStr for ExportFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "obj" => Ok(ExportFormat::Obj),
            "stl" => Ok(ExportFormat::Stl),
            "raw" => Ok(ExportFormat::Raw),
            "glb" => Ok(ExportFormat::Glb),
            _ => Err(format!("Unsupported format: {}", s)),
        }
    }
}

/// Export mesh to the specified format
pub fn export_mesh<W: Write>(mesh: &Mesh, format: ExportFormat, writer: &mut W) -> io::Result<()> {
    match format {
        ExportFormat::Obj => mesh.export_obj(writer),
        ExportFormat::Stl => mesh.export_stl(writer),
        ExportFormat::Raw => mesh.export_raw(writer),
        ExportFormat::Glb => mesh.export_glb(writer),
    }
}

impl Mesh {
    /// Export mesh to OBJ format
    pub fn export_obj<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "# Ocean Generator Mesh")?;
        writeln!(writer, "# Vertices: {}", self.vertices.len())?;
        writeln!(writer, "# Faces: {}", self.faces.len())?;
        
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
        
        // Write faces (OBJ indices are 1-based)
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
    pub fn export_stl<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Write STL header (80 bytes)
        writer.write_all(b"Ocean Generator Mesh                                                        ")?;
        
        // Write number of triangles (4 bytes)
        writer.write_u32::<LittleEndian>(self.faces.len() as u32)?;
        
        // Write triangles
        for face in &self.faces {
            let v1 = self.vertices[face.0].position;
            let v2 = self.vertices[face.1].position;
            let v3 = self.vertices[face.2].position;
            
            let normal = self.vertices[face.0].normal;
            
            // Normal
            writer.write_f32::<LittleEndian>(normal.x)?;
            writer.write_f32::<LittleEndian>(normal.y)?;
            writer.write_f32::<LittleEndian>(normal.z)?;
            
            // Vertices
            writer.write_f32::<LittleEndian>(v1.x)?;
            writer.write_f32::<LittleEndian>(v1.y)?;
            writer.write_f32::<LittleEndian>(v1.z)?;
            
            writer.write_f32::<LittleEndian>(v2.x)?;
            writer.write_f32::<LittleEndian>(v2.y)?;
            writer.write_f32::<LittleEndian>(v2.z)?;
            
            writer.write_f32::<LittleEndian>(v3.x)?;
            writer.write_f32::<LittleEndian>(v3.y)?;
            writer.write_f32::<LittleEndian>(v3.z)?;
            
            // Attribute byte count (2 bytes) - unused
            writer.write_u16::<LittleEndian>(0)?;
        }
        
        Ok(())
    }
    
    /// Export mesh to a simple raw format (vertices and indices)
    pub fn export_raw<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Write header
        writeln!(writer, "Ocean Generator Raw Mesh")?;
        writeln!(writer, "Vertices: {}", self.vertices.len())?;
        writeln!(writer, "Faces: {}", self.faces.len())?;
        
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
    
    /// Export mesh to GLB format (GL Transmission Format Binary) using the gltf crate
    pub fn export_glb<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Create buffers for binary data
        let positions_len = self.vertices.len() * 3 * std::mem::size_of::<f32>();
        let normals_len = self.vertices.len() * 3 * std::mem::size_of::<f32>();
        let texcoords_len = self.vertices.len() * 2 * std::mem::size_of::<f32>();
        let indices_len = self.faces.len() * 3 * std::mem::size_of::<u32>();
        
        let mut buffer_data = Vec::with_capacity(positions_len + normals_len + texcoords_len + indices_len);
        
        // Write positions
        for vertex in &self.vertices {
            buffer_data.write_f32::<LittleEndian>(vertex.position.x)?;
            buffer_data.write_f32::<LittleEndian>(vertex.position.y)?;
            buffer_data.write_f32::<LittleEndian>(vertex.position.z)?;
        }
        
        // Write normals
        for vertex in &self.vertices {
            buffer_data.write_f32::<LittleEndian>(vertex.normal.x)?;
            buffer_data.write_f32::<LittleEndian>(vertex.normal.y)?;
            buffer_data.write_f32::<LittleEndian>(vertex.normal.z)?;
        }
        
        // Write texture coordinates
        for vertex in &self.vertices {
            buffer_data.write_f32::<LittleEndian>(vertex.uv.x)?;
            buffer_data.write_f32::<LittleEndian>(vertex.uv.y)?;
        }
        
        // Write indices
        for face in &self.faces {
            buffer_data.write_u32::<LittleEndian>(face.0 as u32)?;
            buffer_data.write_u32::<LittleEndian>(face.1 as u32)?;
            buffer_data.write_u32::<LittleEndian>(face.2 as u32)?;
        }
        
        // Pad buffer data to multiple of 4 bytes as required by glTF spec
        while buffer_data.len() % 4 != 0 {
            buffer_data.push(0);
        }
        
        // Create simplified GLB export without validation
        let json_chunk = json_export_glb(&self.vertices, &self.faces, buffer_data.len())?;
        
        // GLB header
        writer.write_u32::<LittleEndian>(0x46546C67)?; // Magic: "glTF"
        writer.write_u32::<LittleEndian>(2)?; // Version
        
        let total_length = 12 + 8 + json_chunk.len() + 8 + buffer_data.len();
        writer.write_u32::<LittleEndian>(total_length as u32)?; // Total length
        
        // JSON chunk
        writer.write_u32::<LittleEndian>(json_chunk.len() as u32)?;
        writer.write_u32::<LittleEndian>(0x4E4F534A)?; // Chunk type: "JSON"
        writer.write_all(&json_chunk)?;
        
        // BIN chunk
        writer.write_u32::<LittleEndian>(buffer_data.len() as u32)?;
        writer.write_u32::<LittleEndian>(0x004E4942)?; // Chunk type: "BIN\0"
        writer.write_all(&buffer_data)?;
        
        Ok(())
    }
}

// Helper function to generate GLB JSON without using validation
fn json_export_glb(vertices: &[crate::mesh::Vertex], faces: &[crate::mesh::Face], buffer_length: usize) -> io::Result<Vec<u8>> {
    let vertex_count = vertices.len();
    let positions_byte_length = vertex_count * 3 * 4; // 3 floats per vertex, 4 bytes per float
    let normals_byte_length = vertex_count * 3 * 4;
    let texcoords_byte_length = vertex_count * 2 * 4; // 2 floats per texcoord
    let indices_byte_length = faces.len() * 3 * 4;  // 3 indices per face, 4 bytes per index
    
    let json = serde_json::json!({
        "asset": {
            "version": "2.0",
            "generator": "Ocean Generator Library"
        },
        "scene": 0,
        "scenes": [
            {
                "nodes": [0]
            }
        ],
        "nodes": [
            {
                "mesh": 0,
                "name": "OceanMesh"
            }
        ],
        "meshes": [
            {
                "primitives": [
                    {
                        "attributes": {
                            "POSITION": 0,
                            "NORMAL": 1,
                            "TEXCOORD_0": 2
                        },
                        "indices": 3,
                        "mode": 4 // TRIANGLES
                    }
                ],
                "name": "OceanMesh"
            }
        ],
        "accessors": [
            {
                "bufferView": 0,
                "componentType": 5126, // FLOAT
                "count": vertex_count,
                "type": "VEC3",
                "min": [-1.0, -1.0, -1.0], // Simplified bounding box
                "max": [1.0, 1.0, 1.0]
            },
            {
                "bufferView": 1,
                "componentType": 5126, // FLOAT
                "count": vertex_count,
                "type": "VEC3"
            },
            {
                "bufferView": 2,
                "componentType": 5126, // FLOAT
                "count": vertex_count,
                "type": "VEC2"
            },
            {
                "bufferView": 3,
                "componentType": 5125, // UNSIGNED_INT
                "count": faces.len() * 3,
                "type": "SCALAR"
            }
        ],
        "bufferViews": [
            {
                "buffer": 0,
                "byteLength": positions_byte_length,
                "byteOffset": 0,
                "target": 34962 // ARRAY_BUFFER
            },
            {
                "buffer": 0,
                "byteLength": normals_byte_length,
                "byteOffset": positions_byte_length,
                "target": 34962 // ARRAY_BUFFER
            },
            {
                "buffer": 0,
                "byteLength": texcoords_byte_length,
                "byteOffset": positions_byte_length + normals_byte_length,
                "target": 34962 // ARRAY_BUFFER
            },
            {
                "buffer": 0,
                "byteLength": indices_byte_length,
                "byteOffset": positions_byte_length + normals_byte_length + texcoords_byte_length,
                "target": 34963 // ELEMENT_ARRAY_BUFFER
            }
        ],
        "buffers": [
            {
                "byteLength": buffer_length
            }
        ]
    });
    
    let mut json_string = serde_json::to_string(&json)?;
    
    // Ensure the JSON data is padded to a multiple of 4 bytes
    while json_string.len() % 4 != 0 {
        json_string.push(' ');
    }
    
    Ok(json_string.into_bytes())
}
