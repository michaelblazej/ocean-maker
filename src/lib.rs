// Export main modules
mod mesh;
pub mod wave;
mod export;

// Re-export everything for public use
pub use mesh::{Mesh, Vertex, Face};
pub use wave::WaveParams;
pub use export::{export_mesh, export_mesh_tiled, ExportFormat};

pub mod prelude {
    pub use crate::mesh::{Mesh, Vertex, Face};
    pub use crate::wave::WaveParams;
    pub use crate::export::{export_mesh, export_mesh_tiled, ExportFormat};
}
