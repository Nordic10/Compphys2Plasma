#include "pic_main.h"

int main()
{
  //--------------------Simulation Constants--------------------//
  // Particle Count
  np = 10000;

  // Electron Density
  float ne = 1e15f;
  
  // Grid Number
  nx = 128;
  ny = 128;
  
  // Grid Spacing
  dx = 1e-3f;
  dy = 1e-3f;
  dz = 1e-3f;
  
  // Time Step
  dt = 1e-12f;
  
  // Step Count
  nt = 10000;
  
  // Material Constants
  c = 2.9979e8f;
  e0 = 8.85e-12f;
  mu0 = 1.26e-6f;
  
  // Particle Constants
  q = 1.6e-19f;
  m = 9.1e-31f;

  // Particle Weight
  w = q * ne * (nx * ny) / np;
  m = w * (dx * dy) * m / q; // Mass per Macroparticle
  q = w * (dx * dy); // Charge per Macroparticle
  
  
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
      float xcoord = x-float(nx/2);
      float ycoord = y-float(ny/2);
      float Emag = 1e7 * (50 - sqrt(xcoord*xcoord + ycoord*ycoord));
      float Ex = ((x-nx/2)*Emag);
      float Ey = ((y-ny/2)*Emag);
      //std::cout << "Ex: " << Ex << std::endl;
      //std::cout << "Ey: " << Ey << std::endl;
	    ext_field->E[nx*y+x] = {Ex, Ey, 0};
	    ext_field->B[nx*y+x] = {0, 0, 0};
      // Magnetic field magitude - 1e20
	  }
    }
      
  particle* particles = new particle[np];
  grid_struct* grid = new grid_struct;

  distribute_particles(particles, grid);
  initialize_fields(grid);
  print_output("outfile.txt", particles, false);
  for (int t = 0; t < nt; ++t)
    {
      if (t%100 == 0) 
        std::cout << t << std::endl;
      deposit_fields(particles, grid);
      update_boundary(grid);
      update_fields(grid);
      gather_fields(particles, grid, ext_field);
      push_particles(particles);
      //reset_grid(grid);
      if (t%10)
	print_output("outfile.txt", particles, true);
    }
  print_output("outfile.txt", particles, true);
  print_fields("outfields.txt", ext_field);
  
  return 0;
}

