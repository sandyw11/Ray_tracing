//
//  vector.cpp
//  ascii_image
//
//  Created by Alex on 9/18/16.
//  Copyright © 2016 Alex. All rights reserved.
//

#include "vec3d.h"
#include <cassert>

// Constructors
vec3d::vec3d(float value)
{
    x = y = z = value;
}

vec3d::vec3d(float x, float y, float z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

vec3d::vec3d(const vec3d& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

// Inspectors
float vec3d::operator[](size_t index) const
{
    assert(index < 3);
    return data[index];
}

float& vec3d::operator[](size_t index)
{
    assert(index < 3);
    return data[index];
}

// Operators
vec3d& vec3d::operator=(const vec3d& v)
{
    _assign(v);
    return *this;
}

bool vec3d::operator==(const vec3d& v) const
{
    return (x == v.x) && (y == v.y) && (z == v.z);
}

bool vec3d::operator!=(const vec3d& v) const
{
    return (x != v.x) || (y != v.y) || (z != v.z);
}

vec3d vec3d::operator+(const vec3d& v) const
{
    return vec3d(x + v.x, y + v.y, z + v.z);
}

vec3d vec3d::operator-(const vec3d& v) const
{
    return vec3d(x - v.x, y - v.y, z - v.z);
}

vec3d vec3d::operator*(const vec3d& v) const
{
    return vec3d(x * v.x, y * v.y, z * v.z);
}

vec3d vec3d::operator*(float scale) const
{
    return vec3d(x * scale, y * scale, z * scale);
}

vec3d vec3d::operator/(const vec3d& v) const
{
    return vec3d(x / v.x, y / v.y, z / v.z);
}


vec3d vec3d::operator/(float scale) const
{
    return vec3d(x / scale, y / scale, z / scale);
}

vec3d& vec3d::operator+=(const vec3d& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

vec3d& vec3d::operator-=(const vec3d& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

vec3d& vec3d::operator*=(const vec3d& v)
{
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
}

vec3d& vec3d::operator*=(float scale)
{
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;
}

vec3d& vec3d::operator/=(const vec3d& v)
{
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
}

vec3d& vec3d::operator/=(float scale)
{
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
}

// Modifiers
float vec3d::dot(const vec3d& v) const
{
    return x * v.x + y * v.y + z * v.z;
}

float vec3d::squared_length(void) const
{
    return dot(*this);
}

float vec3d::length(void) const
{
    return sqrt(squared_length());
}

float vec3d::squared_distance(const vec3d& v) const
{
    return (*this - v).squared_length();
}

float vec3d::distance(const vec3d& v) const
{
    return sqrt(squared_distance(v));
}

vec3d vec3d::normalize(void) const
{
    return vec3d(x / length(), y / length(), z / length());
}

vec3d vec3d::cross(const vec3d& v) const
{
    return vec3d(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

// Private Methods
void vec3d::_assign(const vec3d& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

void vec3d::_swap(vec3d& v)
{
    std::swap(x, v.x);
    std::swap(y, v.y);
    std::swap(z, v.z);
}
