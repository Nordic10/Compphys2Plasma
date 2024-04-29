#pragma once
#include "pic_math.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// Simulation Constants
int nx;
int ny;
int nt;
int np;
float dx;
float dy;
float dz;
float dt;
float c;
float e0;
float mu0;
float q;
float m;
float w;


void set_grid(grid_struct* grid)
{
  grid->rho = new float[nx * ny]();
  grid->phi = new float[nx * ny]();
  grid->j = new float3[nx * ny]();
  grid->E = new float3[nx * ny]();
  grid->B = new float3[nx * ny]();
}

void distribute_particles(particle* p, grid_struct* grid)
{
  int x, y, idx;
  float x_, y_;
  set_grid(grid);
  for (int i = 0; i < np; ++i)
    {
      p[i].r.x = random_float(i, true) * dx * nx + dx;
      p[i].r.y = random_float(i, false) * dy * ny + dy;
      x = (int)roundf(p[i].r.x / dx);
      y = (int)roundf(p[i].r.y / dy);
      x_ = abs((p[i].r.x / dx - x) + 0.5f);
      y_ = abs((p[i].r.y / dx - y) + 0.5f);
      idx = nx * y + x;
      grid->rho[idx] += w * y_ * x_;
      grid->rho[idx - 1] += w * y_ * (1 - x_);
      grid->rho[idx - nx] += w * (1 - y_) * x_;
      grid->rho[idx - nx - 1] += w * (1 - y_) * (1 - x_);
    }
}

void initialize_fields(grid_struct* grid)
{
  int idx = 0;
  float tol = 0.0001f;
  float temp = 0;
  bool loop = true;
  float dx2 = dx * dx;
  float dy2 = dy * dy;
  while (loop)
    {
      idx = nx;
      for (int y = 1; y < ny - 1; y++)
	{
	  for (int x = 1; x < nx - 1; x++)
	    {
	      ++idx;
	      if (loop)
		{
		  temp = grid->phi[idx];
		  grid->phi[idx] = 0.5f * (dy2 * (grid->phi[idx + 1] + grid->phi[idx - 1]) + dx2 * (grid->phi[idx + nx] + grid->phi[idx - nx]) + dy2 * dx2 * grid->rho[idx] / e0) / (dy2 + dx2);
		  temp = abs(temp - grid->phi[idx]) / temp;
		  if (temp < tol)
		    loop = false;
		}
	      else
		grid->phi[idx] = 0.5f * (dy2 * (grid->phi[idx + 1] + grid->phi[idx - 1]) + dx2 * (grid->phi[idx + nx] + grid->phi[idx - nx]) + dy2 * dx2 * grid->rho[idx] / e0) / (dy2 + dx2);
	    }
	  idx += 2;
	}
    }
  idx = nx;
  for (int y = 1; y < ny; y++)
    {
      for (int x = 1; x < nx; x++)
	{
	  ++idx;
	  grid->E[idx].x = (grid->phi[idx - 1] - grid->phi[idx]) / dx;
	  grid->E[idx].y = (grid->phi[idx -nx] - grid->phi[idx]) / dy;
	}
      ++idx;
    }
}

void deposit_fields(particle* list, grid_struct* grid)
{
  particle* p;
  float x_, y_, _x, _y;
  int x, y, idx;
  for (int i = 0; i < np; i++)
    {
      p = &list[i];
      x = (int)roundf(p->r.x / dx);
      y = (int)roundf(p->r.y / dy);
      x_ = abs((p->r.x / dx - x) + 0.5f);
      y_ = abs((p->r.y / dx - y) + 0.5f);
      p->r += dt * p->v;
      _x = abs((p->r.x / dx - x) + 0.5f);
      _y = abs((p->r.y / dx - y) + 0.5f);
      idx = nx * y + x;
      grid->j[idx + nx].x += p->v.x * w * (y_ + _y) / 2;
      grid->j[idx + 1].y += p->v.y * w * (x_ + _x) / 2;
      grid->j[idx].x += p->v.x * w * (2 - y_ - _y) / 2;
      grid->j[idx].y += p->v.y * w * (2 - x_ - _x) / 2;
    }
}

void update_boundary(grid_struct* grid)
{
  int i = nx - 1;
  for (int y = 0; y < ny - 1; y++)
    {
      grid->B[i].z = grid->B[i - 1].z;
      i += nx;
    }
  i = nx * (ny - 1);
  for (int x = 0; x < nx; x++)
    {
      grid->B[i].z = grid->B[i - nx].z;
      ++i;
    }
  
  i = 0;
  for (int y = 0; y < ny - 1; y++)
    {
      grid->E[i].x = grid->E[i + 1].x;
      grid->E[i].y = grid->E[i + 1].y;
      i += nx;
    }
  i = 0;
  for (int x = 0; x < nx; x++)
    {
      grid->E[i].x = grid->E[i + nx].x;
      grid->E[i].y = grid->E[i + nx].y;
      ++i;
    }
}

void update_fields(grid_struct* grid)
{
  // B[n-1/2] --> B[n]
  int i = 0;
  for (int y = 0; y < ny - 1; y++)
    {
      for (int x = 0; x < nx - 1; x++)
	{
	  //grid->B[i].x += dt/2 * ((grid->E[i+nxy].y - grid->E[i].y) / dz - (grid->E[i +nx].z - grid->E[i].z) / dy);
	  //grid->B[i].y += dt/2 * ((grid->E[i + 1].z - grid->E[i].z) / dx - (grid->E[i+nxy].x - grid->E[i].x) / dz);
	  grid->B[i].z += dt/2 * ((grid->E[i +nx].x - grid->E[i].x) / dy - (grid->E[i + 1].y - grid->E[i].y) / dx);
	  ++i;
	}
      ++i;
    }
  // E[n-1/2] --> E[n]
  i = nx;
  for (int y = 1; y < ny; y++)
    {
      for (int x = 1; x < nx; x++)
	{
	  ++i;
	  grid->E[i].x -= dt/2 * ((/*(grid->B[i].y - grid->B[i-nxy].y) / dz*/ - (grid->B[i].z - grid->B[i -nx].z) / dy) / mu0 + grid->j[i].x) / e0;
	  grid->E[i].y -= dt/2 * (((grid->B[i].z - grid->B[i - 1].z) / dx /*- (grid->B[i].x - grid->B[i-nxy].x) / dz*/) / mu0 + grid->j[i].y) / e0;
	  //grid->E[i].z -= dt/2 * (((grid->B[i].x - grid->B[i -nx].x) / dy - (grid->B[i].y - grid->B[i - 1].y) / dx) / mu0 + grid->j[i].z) / e0;
	}
      ++i;
    }
  // E[n] --> E[n+1/2]
  i = nx;
  for (int y = 1; y < ny; y++)
    {
      for (int x = 1; x < nx; x++)
	{
	  ++i;
	  grid->E[i].x -= dt/2 * ((/*(grid->B[i].y - grid->B[i-nxy].y) / dz*/ - (grid->B[i].z - grid->B[i -nx].z) / dy) / mu0 + grid->j[i].x) / e0;
	  grid->E[i].y -= dt/2 * (((grid->B[i].z - grid->B[i - 1].z) / dx /*- (grid->B[i].x - grid->B[i-nxy].x) / dz*/) / mu0 + grid->j[i].y) / e0;
	  //grid->E[i].z -= dt/2 * (((grid->B[i].x - grid->B[i -nx].x) / dy - (grid->B[i].y - grid->B[i - 1].y) / dx) / mu0 + grid->j[i].z) / e0;
	}
      ++i;
    }
  // B[n] --> B[n+1/2]
  i = 0;
  for (int y = 0; y < ny - 1; y++)
    {
      for (int x = 0; x < nx - 1; x++)
	{
	  //grid->B[i].x += dt/2 * ((grid->E[i+nxy].y - grid->E[i].y) / dz - (grid->E[i +nx].z - grid->E[i].z) / dy);
	  //grid->B[i].y += dt/2 * ((grid->E[i - 1].z - grid->E[i].z) / dx - (grid->E[i+nxy].x - grid->E[i].x) / dz);
	  grid->B[i].z += dt/2 * ((grid->E[i +nx].x - grid->E[i].x) / dy - (grid->E[i + 1].y - grid->E[i].y) / dx);
	  ++i;
	}
      ++i;
    }
}

void gather_fields(particle* list, grid_struct* grid, external_field* ext)
{
  float x_, y_;
  int x, y, idx;
  float weight;
  particle* p;
  for (int i = 0; i < np; i++)
    {
      p = &list[i];
      x = (int)roundf(p->r.x / dx);
      y = (int)roundf(p->r.y / dy);
      x_ = abs((p->r.x / dx - x) + 0.5f);
      y_ = abs((p->r.y / dx - y) + 0.5f);
      idx = nx * y + x;
      
      if (x == 0 || y == 0)
	{
	  if (x == 0 && y == 0)
	    {
	      p->E.x = (grid->E[idx].x + ext->E[idx].x) * y_;
	      p->E.y = (grid->E[idx].y + ext->E[idx].y) * x_;
	    }
	  else if (x == 0)
	    {
	      p->E.x = (grid->E[idx].x + ext->E[idx].x) * y_ + (grid->E[idx-1].x + ext->E[idx-1].x) * (1 - y_);
	      p->E.y = (grid->E[idx].y + ext->E[idx].y) * x_;
	    }
	  else
	    {
	      p->E.x = (grid->E[idx].x + ext->E[idx].x)* y_;
	      p->E.y = (grid->E[idx].y + ext->E[idx].y) * x_ + (grid->E[idx-nx].y + ext->E[idx-nx].y) * (1 - x_);
	    }
	}
      else
	{
	  p->E.x = (grid->E[idx].x + ext->E[idx].x) * y_ + (grid->E[idx-1].x + ext->E[idx-1].x) * (1 - y_);
	  p->E.y = (grid->E[idx].y + ext->E[idx].y) * x_ + (grid->E[idx-nx].y + ext->E[idx-nx].y) * (1 - x_);
	}
      
      if (x != 0 && y != 0)
	p->B.z = (grid->B[idx].z + ext->B[idx].z) * x_ * y_ + (grid->B[idx-1].z + ext->B[idx-1].z) * (1 - x_) * y_ + (grid->B[idx-nx].z + ext->B[idx-nx].z) * x_ * (1 - y_) + (grid->B[idx-nx-1].z + ext->B[idx-nx-1].z) * (1 - x_) * (1 - y_);
      else if (x == 0 && y != 0)
	p->B.z = (grid->B[idx].z + ext->B[idx].z) * x_ * y_ + (grid->B[idx-nx].z + ext->B[idx-nx].z) * x_ * (1 - y_);
      else if (y == 0 && x != 0)
	p->B.z = (grid->B[idx].z + ext->B[idx].z) * x_ * y_ + (grid->B[idx-1].z + ext->B[idx-1].z) * (1 - x_) * y_;
      else
	p->B.z = (grid->B[idx].z + ext->B[idx].z) * x_ * y_;
    }
}

void push_particles(particle* list)
{
  float3 p0, p1, p_, t, s;
  particle* p; float a, b;
  for (int i = 0; i < np; i++)
    {
      // Standard Boris Pusher
      p = &list[i];
      a = q * dt / 2;
      // If speeds get too close to c this value will become nan
      b = 1 - mag(p->v / c);
      if (b <= 0) 
        b=1e-7;
      p0 = m * p->v / sqrtf(b) + a * p->E;
      t = a * p->B;
      s = 2 * t / (1 + mag(t));
      p1 = p0 + cross(p0 + cross(p0, t), s);
      p_ = p1 + a * p->E;
      p->v = p_ / (m * sqrtf(1 + mag(p_ / (m * c))));
      
      if (p->r.x + dt * p->v.x > nx * dx)
	{
	  p->v.x *= -1;
	  p->r.x = 2 * nx * dx - p->r.x;
	}
      else if (p->r.x + dt * p->v.x < 0)
	{
	  p->v.x *= -1;
	  p->r.x *= -1;
	}
      if (p->r.y + dt * p->v.y > ny * dy)
	{
	  p->v.y *= -1;
	  p->r.y = 2 * ny * dy - p->r.y;
	}
      else if (p->r.y + dt * p->v.y < 0)
	{
	  p->v.y *= -1;
	  p->r.y *= -1;
	}
    }
}

void reset_grid(grid_struct* grid)
{
  for (int i = 0; i < nx * ny; i++)
    grid->j[i] = { 0, 0, 0 };
}

void print_output(const char* filename, particle* p_list,  bool app)
{
  std::ofstream outfile;
  if (app)
    outfile.open(filename, std::ios_base::app);
  else
    outfile.open(filename);
  particle p;
  for (int i = 0; i < np; i++)
    {
      p = p_list[i];
      outfile << p.r.x << ',' << p.r.y << ',' << p.r.z << std::endl;
    } 
}

void print_fields(const char* filename, external_field* ext_field) {
  std::ofstream outfile;
  outfile.open(filename);
  float3 E;
  float3 B;

  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      E = ext_field->E[nx*y+x];
      B = ext_field->B[nx*y+x];

      outfile << x << "," << y << "," << E.x << "," << E.y << "," << E.z << "," << B.x << "," << B.y << "," << B.z << std::endl;
    }
  }
}

float3* read_Efields(const char* filename) {
  std::ifstream infile(filename);
  std::string line;

    float3* E = new float3[nx*ny]();

  while (std::getline(infile, line)) {
    int x, y;
    float dx, dy;

    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> row;
    while (std::getline(ss, cell, ',')) {
      row.push_back(cell);
    }
    x = std::stoi(row[0]);
    y = std::stoi(row[1]);
    dx = 1e7 * std::stof(row[2]);
    dy = 1e7 * std::stof(row[3]);

    E[nx*y+x] = {dx,dy, 0};
  }

  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      if (i==0 || i==nx-1 || j==0 || j==ny-1) {
        E[nx*j+i] = {0,0,0}; //gradient doesn't cover electric field on boundaries, so they are set to 0
      }
    }
  }

  return E;
}

