#pragma once
#include "pic_structs.h"

float noise(float x, float y, int seed)
{
	const unsigned w = 8 * sizeof(unsigned);
	const unsigned s = w / 2;
	unsigned a = x, b = y, c = seed;
	a *= 3284157443; b ^= a << s | a >> w - s;
	b *= 1911520717; a ^= b << s | b >> w - s;
	c *= 4628058208; a ^= c << s | c >> w - s;
	a *= 2048419325;
	float random = a * (0.5 / ~(~0u >> 1));
	return random;
}

float random_float(int seed, bool xy)
{
	return noise(3*(seed+xy), ~(2-xy*1137)/seed, seed);
}

float mag(float3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

float dot(float3 v1, float3 v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

float3 div(float3 v1, float3 v2)
{
	return { v1.x / v2.x, v1.y / v2.y, v1.z / v2.z };
}

float3 cross(float3 v1, float3 v2)
{
	return { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
}

float3& operator+=(float3& lhs, const float3& rhs)
{
	lhs.x += rhs.x;
	lhs.y += rhs.y;
	lhs.z += rhs.z;
	return lhs;
}

float3 operator+(const float3& lhs, const float3& rhs)
{
	float3 result;
	result.x = lhs.x + rhs.x;
	result.y = lhs.y + rhs.y;
	result.z = lhs.z + rhs.z;
	return result;
}

float3 operator-(const float3& lhs, const float3& rhs)
{
	float3 result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.z = lhs.z - rhs.z;
	return result;
}

float3 operator*(const float3& lhs, const float& rhs)
{
	float3 result;
	result.x = lhs.x * rhs;
	result.y = lhs.y * rhs;
	result.z = lhs.z * rhs;
	return result;
}
float3 operator*(const float& lhs, const float3& rhs)
{
	float3 result;
	result.x = lhs * rhs.x;
	result.y = lhs * rhs.y;
	result.z = lhs * rhs.z;
	return result;
}

float3 operator/(const float3& lhs, const float& rhs)
{
	float3 result;
	result.x = lhs.x / rhs;
	result.y = lhs.y / rhs;
	result.z = lhs.z / rhs;
	return result;
}
