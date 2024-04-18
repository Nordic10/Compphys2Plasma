#pragma once
#include "pic_math.h"
#include <math.h>
#include <fstream>

void set_grid(grid_struct* grid, int nx, int ny)
{
  grid->rho = new float[nx * ny]();
  grid->phi = new float[nx * ny]();
  grid->j = new float3[nx * ny]();
  grid->E = new float3[nx * ny]();
  grid->B = new float3[nx * ny]();
}

void distribute_particles(particle* p, grid_struct* grid, int np, int nx, int ny, float dx, float dy, float q)
{
  int x, y, idx;
  float x_, y_, weight;
  set_grid(grid, nx, ny);
  for (int i = 0; i < np; i++)
    {
      p[i].r.x = random_float(i, true) * dx;
      p[i].r.y = random_float(i, false)* dy;
      x = (int)roundf(p[i].r.x / dx);
      y = (int)roundf(p[i].r.y / dy);
      x_ = abs((p[i].r.x / dx - x) * 2 + 1);
      y_ = abs((p[i].r.y / dx - y) * 2 + 1);
      idx = nx * y + x;
      weight = q;
      grid->rho[idx] += weight * y_ * x_;
      grid->rho[idx - 1] += weight * y_ * (1 - x_);
      grid->rho[idx - nx] += dx * weight * (1 - y_) * x_;
      grid->rho[idx - nx - 1] += dy * weight * (1 - y_) * (1 - x_);
    }
}

void initialize_fields(grid_struct* grid, int nx, int ny, float dx, float dy, float e0)
{
  int idx = nx;
  float tol = 0.1f;
  float temp = 0;
  bool loop = true;
  while (loop)
    {
      for (int y = 1; y < ny - 1; y++)
	{
	  for (int x = 1; x < nx - 1; x++)
	    {
	      ++idx;
	      if (loop)
		{
		  temp = grid->phi[idx];
		  grid->phi[idx] = 0.25f * (grid->phi[idx + 1] + grid->phi[idx - 1] + grid->phi[idx + nx] + grid->phi[idx - nx] + grid->rho[idx] / e0);
		  temp = abs(temp - grid->phi[idx]) / temp;
		  if (temp < tol)
		    loop = false;
		}
	      else
		grid->phi[idx] = 0.25f * (grid->phi[idx + 1] + grid->phi[idx - 1] + grid->phi[idx + nx] + grid->phi[idx - nx] + grid->rho[idx] / e0);
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
    }
}

void field_deposition(particle* list, grid_struct* grid, int np, int nx, int ny, float dx, float dy, float dt, float q)
{
  particle* p;
  float weight;
  float x_, y_;
  int x, y, idx;
  for (int i = 0; i < np; i++)
    {
      p = &list[i];
      x = (int)roundf(p->r.x / dx);
      y = (int)roundf(p->r.y / dy);
      x_ = abs((p->r.x / dx - x) * 2 + 1);
      y_ = abs((p->r.y / dx - y) * 2 + 1);
      idx = nx * y + x;
      weight = q;
      grid->j[idx].x += p->v.x * dt * weight * x_;
      grid->j[idx].y += p->v.y * dt * weight * y_;
      grid->j[idx - 1].x += p->v.x * dt * weight * (1 - x_);
      grid->j[idx -nx].y += p->v.y * dt * weight * (1 - y_);
    }
}

void update_fields(grid_struct* grid, int nx, int ny, float dx, float dy, float dz, float dt, float e0, float mu0)
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
	  grid->E[i].x -= (dt/2 * (/*(grid->B[i].y - grid->B[i-nxy].y) / dz*/ - (grid->B[i].z - grid->B[i -nx].z) / dy) / mu0 + grid->j[i].x) / e0;
	  grid->E[i].y -= (dt/2 * ((grid->B[i].z - grid->B[i - 1].z) / dx /*- (grid->B[i].x - grid->B[i-nxy].x) / dz*/) / mu0 + grid->j[i].y) / e0;
	  //grid->E[i].z -= (dt/2 * ((grid->B[i].x - grid->B[i -nx].x) / dy - (grid->B[i].y - grid->B[i - 1].y) / dx) / mu0 + grid->j[i].z) / e0;
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
	  grid->E[i].x -= (dt/2 * (/*(grid->B[i].y - grid->B[i-nxy].y) / dz*/ - (grid->B[i].z - grid->B[i -nx].z) / dy) / mu0 + grid->j[i].x) / e0;
	  grid->E[i].y -= (dt/2 * ((grid->B[i].z - grid->B[i - 1].z) / dx /*- (grid->B[i].x - grid->B[i-nxy].x) / dz*/) / mu0 + grid->j[i].y) / e0;
	  //grid->E[i].z -= (dt/2 * ((grid->B[i].x - grid->B[i -nx].x) / dy - (grid->B[i].y - grid->B[i - 1].y) / dx) / mu0 + grid->j[i].z) / e0;
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

void field_gathering(particle* list, grid_struct* grid, int np, int nx, int ny, float dx, float dy, float q)
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
      x_ = abs((p->r.x / dx - x) * 2 + 1);
      y_ = abs((p->r.y / dx - y) * 2 + 1);
      idx = nx * y + x;
      weight = q;
      
      if (x == 0 || y == 0)
	{
	  if (x == 0 && y == 0)
	    {
	      p->E.x = weight * grid->E[idx].x * x_;
	      p->E.y = weight * grid->E[idx].y * y_;
	    }
	  else if (x == 0)
	    {
	      p->E.x = weight * (grid->E[idx].x * x_ + grid->E[idx - 1].x * (1 - x_));
	      p->E.y = weight * grid->E[idx].y * y_;
	    }
	  else
	    {
	      p->E.x = weight * grid->E[idx].x * x_;
	      p->E.y = weight * (grid->E[idx].y * y_ + grid->E[idx -nx].y * (1 - y_));
	    }
	}
      else
	{
	  p->E.x = weight * (grid->E[idx].x * x_ + grid->E[idx - 1].x * (1 - x_));
	  p->E.y = weight * (grid->E[idx].y * y_ + grid->E[idx -nx].y * (1 - y_));
	}
      
      if (x != 0 && y != 0)
	p->B.z = weight * (grid->B[idx].z * x_ * y_ + grid->B[idx - 1].x * x_ * (1 - y_) + grid->B[idx - nx].z * (1 - x_) * y_ + grid->B[idx - 1 - nx].x * (1 - x_) * (1 - y_));
      else if (x == 0 && y != 0)
	p->B.z = weight * (grid->B[idx].z * x_ * y_ + grid->B[idx - nx].z * (1 - x_) * y_);
      else if (y == 0 && x != 0)
	p->B.z = weight * (grid->B[idx].z * x_ * y_ + grid->B[idx - 1].x * x_ * (1 - y_));
      else
	p->B.z = weight * (grid->B[idx].z * x_ * y_);
    }
}

void push_particles(particle* list, int np, int nx, int ny, float dt, float q, float m, float c)
{
  float3 p0, p1, p_, t, s;
  particle* p; float a;
  for (int i = 0; i < np; i++)
    {
      p = &list[i];
      a = q * dt / 2;
      p0 = m * p->v / sqrtf(1 - mag(p->v / c)) + a * p->E;
      t = a * p->B;
      s = 2 * t / (1 + mag(t));
      p1 = p0 + cross(p0 + cross(p0, t), s);
      p_ = p1 + a * p->E;
      p->v = p_ / (m * sqrtf(1 + mag(p_ / (m * c))));
      
      if (p->r.x + dt/2 * p->v.x > nx)
	{
	  p->v.x *= -1;
	  p->r.x = 2 * nx - p->r.x;
	}
      else if (p->r.x + dt/2 * p->v.x < 0)
	{
	  p->v.x *= -1;
	  p->r.x *= -1;
	}
      if (p->r.y + dt/2 * p->v.y > ny)
	{
	  p->v.y *= -1;
	  p->r.y = 2 * ny - p->r.y;
	}
      else if (p->r.y + dt/2 * p->v.y < 0)
	{
	  p->v.y *= -1;
	  p->r.y *= -1;
	}
      
      p->r += dt/2 * p->v;
      
      if (p->r.x + dt / 2 * p->v.x > nx)
	{
	  p->v.x *= -1;
	  p->r.x = 2 * nx - p->r.x;
	}
      else if (p->r.x + dt / 2 * p->v.x < 0)
	{
	  p->v.x *= -1;
	  p->r.x *= -1;
	}
      if (p->r.y + dt / 2 * p->v.y > ny)
	{
	  p->v.y *= -1;
	  p->r.y = 2 * ny - p->r.y;
	}
      else if (p->r.y + dt / 2 * p->v.y < 0)
	{
	  p->v.y *= -1;
	  p->r.y *= -1;
	}
      
      p->r += dt / 2 * p->v;
    }
}

void reset_grid(grid_struct* grid, int nx, int ny)
{
  for (int i = 0; i < nx * ny; i++)
    {
      grid->rho[i] = 0;
      grid->j[i] = { 0, 0, 0 };
    }
}

void print_output(const char* filename, particle* p_list, int n, bool app)
{
  std::ofstream outfile;
  if (app)
    outfile.open(filename, std::ios_base::app);
  else
    outfile.open(filename);
  particle p;
  for (int i = 0; i < n; i++)
    {
      p = p_list[i];
      outfile << p.r.x << ',' << p.r.y << ',' << p.r.z << ',' << std::endl;
    } 
}

