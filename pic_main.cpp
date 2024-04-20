#include "pic_main.h"

int main()
{
  //--------------------Simulation Constants--------------------//
  // Particle Count
  np = 10000;

  // Electron Density
  float ne = 1 * powf(10, 6);
  
  // Grid Number
  nx = 128;
  ny = 128;
  
  // Grid Spacing
  dx = 1 * powf(10, -6);
  dy = 1 * powf(10, -6);
  dz = 1 * powf(10, -6);
  
  // Time Step
  dt = 1 * powf(10, -9);
  
  // Step Count
  nt = 100;
  
  // Material Constants
  c = 2.9979f * powf(10, -8);
  e0 = 8.85f * powf(10, -12);
  mu0 = 1.26f * powf(10, -6);
  
  // Particle Constants
  q = 1.6f * powf(10, -19) * (ne * dx * dy / 4);
  m = 9.1f * powf(10, -31) * (ne * dx * dy / 4);

  //--------------------External Fields and Boundary Conditions--------------------//
  //Create External Fields
  external_field* ext_field = new external_field;
  if (ext_field != nullptr)
    {
      ext_field->E = new float3[nx*ny]();
      ext_field->B = new float3[nx*ny]();
      for (int x = 0; x < nx; x++)
	for (int y = 0; y < ny; y++)
	  {
	    //ext_fields->E[y*nx+x] = {f(x), f(y), f(z)};
	    //ext_fields->B[y*nx+x] = {f(x), f(y), f(z)};
	  }
    }

  // Create Walls
  wall_boundary* walls = nullptr;
  if (walls != nullptr)
    for (int x = 0; x < nx; x++)
      for (int y = 0; y < ny; y++)
	{
	  //walls->wall[y*nx+x] = 1;
	}
  
  particle* particles = new particle[np];
  grid_struct* grid = new grid_struct;
  
  distribute_particles(particles, grid);
  initialize_fields(grid);
  print_output("outfile.txt", particles, false);
  for (int t = 0; t < nt; ++t)
    {
      field_deposition(particles, grid);
      update_fields(grid);
      field_gathering(particles, grid, ext_field);
      push_particles(particles);
      reset_grid(grid);
      if (t % 10 == 0)
	print_output("outfile.txt", particles, true);
    }
  
  return 0;
}

