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

struct laser
{
	float E0;
	float w0;
	float k;
	float omega;
	float phi;
	float theta;
	int edge;
	float r0;
	float t0;
};

struct grid_struct
{
	float* rho;
	float* phi;
	float3* j;
	float3* E;
	float3* B;
};
