//
//  vector.h
//  ascii_image
//
//  Created by Alex on 9/18/16.
//  Copyright © 2016 Alex. All rights reserved.
//

#ifndef vec3d_h
#define vec3d_h

#include <iostream>
#include <cmath>

class vec3d {
public:
    union {
        struct { float x, y, z; };
        float data[3];
    };
    
public:
    // Constructor
    explicit vec3d(float value=0.0f);
    vec3d(float x, float y, float z);
    vec3d(const vec3d& v);
    
    // Desctructor
    ~vec3d(void) {}
    
    // Inspectors
    size_t size(void) const { return 3; }
    
    float  operator[](size_t index) const;
    float& operator[](size_t index);
    
    // Operators
    vec3d& operator=(const vec3d& v);
    
    bool operator==(const vec3d& v) const;
    bool operator!=(const vec3d& v) const;
    
    vec3d operator+(const vec3d& v) const;
    vec3d operator-(const vec3d& v) const;
    vec3d operator*(const vec3d& v) const;
    vec3d operator*(float scale) const;
    vec3d operator/(const vec3d& v) const;
    vec3d operator/(float scale) const;
    
    vec3d& operator+=(const vec3d& v);
    vec3d& operator-=(const vec3d& v);
    vec3d& operator*=(const vec3d& v);
    vec3d& operator*=(float scale);
    vec3d& operator/=(const vec3d& v);
    vec3d& operator/=(float scale);
    
    // Modifiers
    float dot(const vec3d& v) const;
    float squared_length(void) const;
    float length(void) const;
    float squared_distance(const vec3d& v) const;
    float distance(const vec3d& v) const;
    
    vec3d normalize(void) const;
    vec3d cross(const vec3d& v) const;
    
    // Friends
    friend void  swap(vec3d& a, vec3d& b)                { return a._swap(b); }
    friend void  normalize(vec3d& v)                     { v /= v.length(); }
    friend vec3d operator*(float scale, const vec3d& v)  { return (v*scale); }
    
    friend std::ostream& operator<<(std::ostream& s, const vec3d& v)
    {
        s << "[" << v.x << "," << v.y << "," << v.z << "]";
        return s;
    }
    
private:
    // Private Methods
    void _assign(const vec3d& v);
    void _swap(vec3d& v);
};

#endif /* vec3d_h */
