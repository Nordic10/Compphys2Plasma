#include "pic_main.h"

int main()
{
  // Particle Count
  int np = 10000;

  // Electron Density
  float ne = 1 * powf(10, 6);
  
  // Laser Count
  int nl = 1;
  
  // Grid Number
  int nx = 127;
  int ny = 127;
  nx++; ny++;
  
  // Grid Spacing
  float dx = 1 * powf(10, -6);
  float dy = 1 * powf(10, -6);
  float dz = 1 * powf(10, -6);
  
  // Time Step
  float dt = 1 * powf(10, -9);
  
  // Step Count
  int nt = 100;
  
  // Material Constants
  float c = 2.9979f * powf(10, -8);
  float e0 = 8.85f * powf(10, -12);
  float mu0 = 1.26f * powf(10, -6);
  
  // Particle Constants
  float q = 1.6f * powf(10, -19) * (ne * dx * dy / 4);
  float m = 9.1f * powf(10, -31) * (ne * dx * dy / 4);
  
  particle* particles = new particle[np];
  grid_struct* grid = new grid_struct;
  
  distribute_particles(particles, grid, np, nx, ny, dx, dy, q);
  initialize_fields(grid, nx, ny, dx, dy, e0);
  print_output("outfile.txt", particles, np, false);
  for (int t = 0; t < nt; ++t)
    {
      field_deposition(particles, grid, np, nx, ny, dx, dy, dt, q);
      //update_lasers(lasers, grid, nl, nx, ny, dx, dy, dz, t);
      update_fields(grid, nx, ny, dx, dy, dz, dt, e0, mu0);
      field_gathering(particles, grid, np, nx, ny, dx, dy, q);
      push_particles(particles, np, nx, ny, dt, q, m, c);
      reset_grid(grid, nx, ny);
      if (t % 10 == 0)
	print_output("outfile.txt", particles, np, true);
    }
  
  return 0;
}

