use glam::Vec2;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::f32::consts::PI;

/// Parameters for a single wave component
#[derive(Debug, Clone, Copy)]
pub struct WaveParams {
    pub amplitude: f32,
    pub wavelength: f32,
    pub steepness: f32,
    pub direction: Vec2,
    pub frequency: f32,
}

/// Generate random wave parameters
pub fn generate_wave_params(
    count: usize,
    base_amplitude: f32,
    base_wavelength: f32,
    base_steepness: f32,
    base_direction_degrees: f32,
    seed: u64,
) -> Vec<WaveParams> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut waves = Vec::with_capacity(count);
    
    // Convert direction from degrees to radians
    let direction_angle = base_direction_degrees * PI / 180.0;
    
    for i in 0..count {
        // Calculate wavelength - shorter for higher frequency components
        let wavelength_factor = 1.0 / (1.0 + i as f32 * 0.2);
        let wavelength = base_wavelength * wavelength_factor;
        
        // Amplitude decreases for higher frequency components
        let amplitude_factor = f32::exp(-0.3 * i as f32);
        let amplitude = base_amplitude * amplitude_factor;
        
        let direction = Vec2::new(f32::cos(direction_angle), f32::sin(direction_angle));
        
        // Randomize steepness a bit
        let steepness = base_steepness * rng.gen_range(0.8..1.2);
        
        // Calculate angular frequency (2π/T)
        let frequency = 2.0 * PI * f32::sqrt(9.8 / wavelength);
        
        waves.push(WaveParams {
            amplitude,
            wavelength,
            steepness: steepness.min(1.0), // Ensure steepness is at most 1.0
            direction,
            frequency,
        });
    }
    
    waves
}
