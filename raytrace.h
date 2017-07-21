//
//  raytrace.h
//  ascii_image
//
//  Created by Alex on 9/18/16.
//  Copyright © 2016 Alex. All rights reserved.
//

#ifndef raytrace_h
#define raytrace_h

#include <stdio.h>
#include "vec3d.h"
#include <cmath>
#include <vector>

using namespace std;

typedef struct {
    vec3d rorg;
    vec3d rdir;
    float eta;
} RayType;

typedef struct {
    int idx;
    vec3d dcolor;
    vec3d scolor;
    float ka;
    float kd;
    float ks;
    int n;
    float alpha;
    float eta;
} MaterialType;

typedef struct {
    int idx;
    vec3d sorg;
    float r;
    int material;
    int texture;
} SphereType;

typedef struct {
    int idx;
    vec3d lorg;
    int w;
    vec3d lcolor;
} LightType;

typedef struct {
    int idx;
    int v1;
    int v2;
    int v3;
    int vt1;
    int vt2;
    int vt3;
    int vn1;
    int vn2;
    int vn3;
    vec3d vertex1;
    vec3d vertex2;
    vec3d vertex3;
    vec3d vertexture1;
    vec3d vertexture2;
    vec3d vertexture3;
    int material;
    int texture;
    bool smooth;
    vec3d normal;
    vec3d normal1;
    vec3d normal2;
    vec3d normal3;
    float alpha;
    float beta;
    float gamma;
} TriangleType;

typedef struct {
    int idx;
    int width;
    int height;
    vec3d texcolor[100][100];
} TextureType;

typedef struct {
    int object_idx;
    int object_type;
    int subobject_idx;
    int material_idx;
    int texture_idx;
} SceneType;

typedef struct {
    int object_idx;
    float distance;
    int object_type;
    int subobject_idx;
    int material_idx;
    int texture_idx;
} RayPayType;

extern void Get_Ray(const RayType, const vec3d, const int, const int, const int, const int, const float, const bool, RayType &);
extern void Trace_Ray(const RayType, vector<SceneType> &, vector<SphereType> &, vector<TriangleType> &, RayPayType &);
extern vec3d Shade_Ray(const RayType, const RayPayType raypay, vector<SceneType> &, vector<SphereType> &, vector<TriangleType> &, vector<MaterialType> &, vector<TextureType> &, vector<LightType> &, vec3d &, int &);

#endif /* raytrace_h */
