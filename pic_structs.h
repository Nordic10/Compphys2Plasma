#pragma once

struct float3
{
  float x = 0;
  float y = 0;
  float z = 0;
};

struct particle
{
  float3 r;
  float3 v;
  float3 E;
  float3 B;
};

struct external_field
{
  float3* E;
  float3* B;
};

struct wall_boundary
{
  bool* wall;
};

struct grid_struct
{
  float* rho;
  float* phi;
  float3* j;
  float3* E;
  float3* B;
};
