//
//  raytrace.cpp
//  ascii_image
//
//  Created by Alex on 9/18/16.
//  Copyright © 2016 Alex. All rights reserved.
//

#include "raytrace.h"

using namespace std;

bool Sphere_Inter(const RayType, const SphereType, float &);
bool Triangle_Inter(const RayType, TriangleType &, float &);
void Shade_Ray(const vec3d, const SphereType, vector<SphereType> &, const int, vector<LightType> &, const int lnum, const RayType, const MaterialType, const TextureType);
void Sphere_Texture(const SphereType, const vec3d, float &, float &);
void Triangle_Texture(const TriangleType, float &, float &);

// check whether there is intersect, whether update the color of the pixel of output image
void Get_Ray(const RayType ray, const vec3d up, const int imgi, const int imgj, const int width, const int height, const float fovv, const bool parallel, RayType & viewray)
{
    vec3d u, v;
    float aspect;
    vec3d n;
    float d, w, h;
    vec3d ul, ur, ll, lr;
    vec3d dh, dv;
    vec3d vpoint;
    float pi = 4.0f * atan(1.0f);
    float radians;
    
    u = ray.rdir.cross(up).normalize();
    v = u.cross(ray.rdir).normalize();
    aspect = float(width) / float(height);
    n = ray.rdir.normalize();
    d = 10.0f;
    radians = fovv / 180.0f * pi;
    h = tan(radians / 2.0f) * d * 2.0f;
    w = h * aspect;
    
    // calculate the four corners of the viewing window
    if(parallel == true)
    {
        ul = ray.rorg + v * h / 2.0f - u * w / 2.0f;
        ur = ray.rorg + v * h / 2.0f + u * w / 2.0f;
        ll = ray.rorg - v * h / 2.0f - u * w / 2.0f;
        lr = ray.rorg - v * h / 2.0f + u * w / 2.0f;
    }
    else
    {
        ul = ray.rorg + n * d + v * h / 2.0f - u * w / 2.0f;
        ur = ray.rorg + n * d + v * h / 2.0f + u * w / 2.0f;
        ll = ray.rorg + n * d - v * h / 2.0f - u * w / 2.0f;
        lr = ray.rorg + n * d - v * h / 2.0f + u * w / 2.0f;
    }
    
    // calculate the position of the point in the viewing window corresponding to the output image
    dh = (ur - ul) / (width - 1.0f);
    dv = (ll - ul) / (height - 1.0f);
    vpoint = ul + dv * imgi + dh * imgj;
    
    viewray.rdir = (vpoint - ray.rorg).normalize();
    viewray.rorg = ray.rorg;
}

void Trace_Ray(const RayType viewray, vector<SceneType> & scene, vector<SphereType> & sphere, vector<TriangleType> & triangle, RayPayType & raypay)
{
    bool intersect = false;
    float distance = 1000.0f;
    float t = -1.0f;
    vector<SceneType>::iterator sciter;
    vector<SphereType>::iterator spiter = sphere.begin();
    vector<TriangleType>::iterator triter = triangle.begin();
    for(sciter = scene.begin(); sciter != scene.end(); sciter++)
    {
        if((*sciter).object_type == 0)
        {
            intersect = Sphere_Inter(viewray, *(spiter + (*sciter).subobject_idx), t);
        }
        else if((*sciter).object_type == 1)
        {
            intersect = Triangle_Inter(viewray, *(triter + (*sciter).subobject_idx), t);
        }
        //if(intersect == true)
        if(t > 0.0f)
        {
            if(t < distance)
            {
                distance = t;
                raypay.distance = t;
                raypay.object_idx = (*sciter).object_idx;
                raypay.object_type = (*sciter).object_type;
                raypay.subobject_idx = (*sciter).subobject_idx;
                raypay.material_idx = (*sciter).material_idx;
                raypay.texture_idx = (*sciter).texture_idx;
            }
        }
        else
        {
            raypay.object_idx = -1;
        }
    }
}

bool Sphere_Inter(const RayType ray, const SphereType sphere, float & t)
{
    bool intersect = false;
    float disc = -1.0f, t0 = -1.0f, t1 = -1.0f, t2 = -1.0f;
    float bias = 1.0e-06;
    // calculate discriminant to check whether current view ray intersect with shpere
    disc = pow(2.0f * ray.rdir.dot(ray.rorg - sphere.sorg), 2) - 4.0f * ((ray.rorg - sphere.sorg).squared_length() - pow(sphere.r, 2));
    
    // if discriminant >= 0, there is intersect
    if(disc > 0)
    {
        t1 = (- 2.0f * ray.rdir.dot(ray.rorg - sphere.sorg) + sqrt(disc)) / 2.0f;
        t2 = (- 2.0f * ray.rdir.dot(ray.rorg - sphere.sorg) - sqrt(disc)) / 2.0f;
        
        // if t > 0, the intersect is in front of eye
        if(t1 > bias)
        {
            //cout <<"t1=" <<t1 <<endl;
            intersect = true;
            t = t1;
            if(t2 > bias && t2 < t1)
            {
                //cout <<"t2=" <<t2 <<endl;
                t = t2;
            }
        }
    }
    else if(disc == 0)
    {
        t0 = -ray.rdir.dot(ray.rorg - sphere.sorg);
        if(t0 > bias)
        {
            //cout <<"t0=" <<t0 <<endl;
            intersect = true;
            t = t0;
        }
    }
    return intersect;
}

bool Triangle_Inter(const RayType ray, TriangleType & triangle, float & t)
{
    bool intersect = false;
    vec3d p0, p1, p2, x0, xd;
    p0 = triangle.vertex1;
    p1 = triangle.vertex2;
    p2 = triangle.vertex3;
    x0 = ray.rorg;
    xd = ray.rdir;
    vec3d e1, e2, n;
    float D, t0;
    vec3d p;
    vec3d e3, e4;
    float A, a, b, c;
    float alpha, beta, gamma;
    float epsilon = 0.00001f;
    e1 = p1 - p0;
    e2 = p2 - p0;
    n = e1.cross(e2);
    triangle.normal = n;
    if(n.dot(xd) != 0.0f)
    {
        D = - n.dot(p0);
        t0 = - (n.dot(x0) + D) / n.dot(xd);
        p = x0 + t0 * xd;
        e3 = p - p1;
        e4 = p - p2;
        A = 0.5f * (e1.cross(e2)).length();
        a = 0.5f * (e3.cross(e4)).length();
        b = 0.5f * (e4.cross(e2)).length();
        c = 0.5f * (e1.cross(e3)).length();
        alpha = a / A;
        beta = b / A;
        gamma = c / A;
        triangle.alpha = alpha;
        triangle.beta = beta;
        triangle.gamma = gamma;
        if(alpha >= 0.0f && alpha <= 1.0f && beta >= 0.0f && beta <= 1.0f && gamma >= 0.0f && gamma <= 1.0f && alpha + beta + gamma - 1.0f < epsilon)
        {
            intersect = true;
            t = t0;
        }
    }
    return intersect;
}

// generate the color of each pixel of the output image
vec3d Shade_Ray(const RayType viewray, const RayPayType raypay, vector<SceneType> & scene, vector<SphereType> & sphere, vector<TriangleType> & triangle, vector<MaterialType> & material, vector<TextureType> & texture, vector<LightType> & light, vec3d & img, int & depth)
{
    vec3d intersect;
    intersect = viewray.rorg + viewray.rdir * raypay.distance;
    vector<MaterialType>::const_iterator mtiter = material.begin() + raypay.material_idx;
    vec3d Od, Os;
    Od = (*mtiter).dcolor;
    Os = (*mtiter).scolor;
    
    vector<SphereType>::iterator spiter;
    vector<TriangleType>::iterator triter;
    vec3d N;
    if(raypay.object_type == 0)
    {
        spiter = sphere.begin() + raypay.subobject_idx;
        N = (intersect - (*spiter).sorg) / (*spiter).r;
        
        if(raypay.texture_idx > -1)
        {
            float u, v;
            int i, j;
            Sphere_Texture(*spiter, intersect, u, v);
            vector<TextureType>::const_iterator texiter = texture.begin() + raypay.texture_idx;
            i = round(u * (float)((*texiter).width - 1));
            j = round(v * (float)((*texiter).height - 1));
            Od = (*texiter).texcolor[i][j];
        }
    }
    else if(raypay.object_type == 1)
    {
        triter = triangle.begin() + raypay.subobject_idx;
        if((*triter).smooth == true)
        {
            float alpha, beta, gamma;
            vec3d n0, n1, n2;
            alpha = (*triter).alpha;
            beta = (*triter).beta;
            gamma = (*triter).gamma;
            n0 = (*triter).normal1;
            n1 = (*triter).normal2;
            n2 = (*triter).normal3;
            N = (alpha * n0 + beta * n1 + gamma * n2) / (alpha * n0 + beta * n1 + gamma * n2).length();
        }
        else
        {
            N = (*triter).normal;
        }
        
        if(raypay.texture_idx > -1)
        {
            float u, v;
            int i, j;
            Triangle_Texture(*triter, u, v);
            vector<TextureType>::const_iterator texiter = texture.begin() + raypay.texture_idx;
            i = round(u * (float)((*texiter).width - 1));
            j = round(v * (float)((*texiter).height - 1));
            Od = (*texiter).texcolor[i][j];
        }
    }
    
    vec3d L, V, H;
    V = viewray.rdir * (-1.0f);
    float NL = 0.0f, NH = 0.0f;
    vec3d I;
    RayPayType lraypay;
    float bias = 1.0e-06;
    vector<LightType>::const_iterator lgiter;
    for(lgiter = light.begin(); lgiter != light.end(); lgiter++)
    {
        RayType sharay;
        int flag = 1;
        
        if((*lgiter).w == 0)
        {
            L = (*lgiter).lorg.normalize() * (-1.0f);
            H = (L + V).normalize();
            NL = N.dot(L);
            if(NL < 0.0f)
            {
                NL = 0.0f;
            }
            NH = N.dot(H);
            if(NH < 0.0f)
            {
                NH = 0.0f;
            }
            sharay.rorg = intersect;
            sharay.rdir = L;
            
            Trace_Ray(sharay, scene, sphere, triangle, lraypay);
            
            if(lraypay.object_idx > -1)
            {
                if(lraypay.distance > bias)
                {
                    flag = 0;
                }
            }
        }
        else if((*lgiter).w == 1)
        {
            float dis = 0.0f;
            L = ((*lgiter).lorg - intersect).normalize();
            dis = ((*lgiter).lorg - intersect).length();
            H = (L + V).normalize();
            NL = N.dot(L);
            if(NL < 0.0f)
            {
                NL = 0.0f;
            }
            NH = N.dot(H);
            if(NH < 0.0f)
            {
                NH = 0.0f;
            }
            sharay.rorg = intersect;
            sharay.rdir = L;
            
            Trace_Ray(sharay, scene, sphere, triangle, lraypay);
            
            if(lraypay.object_idx > -1)
            {
                if(lraypay.distance < dis && lraypay.distance > bias)
                {
                    flag = 0;
                }
            }
        }
        //I += flag * (*lgiter).lcolor * ((*mtiter).kd * Od * NL + (*mtiter).ks * Os * pow(NH, (*mtiter).n));
        I += flag * (1 - (*mtiter).alpha) * (*lgiter).lcolor * ((*mtiter).kd * Od * NL + (*mtiter).ks * Os * pow(NH, (*mtiter).n));
    }
    I += (*mtiter).ka * Od;
    
    //cout <<"org depth=" <<depth <<endl;
    if(depth < 5)
    {
        float Fr, F0, etai, etat;
        vec3d reflecolor(0.0f);
        int refledepth;
        refledepth = depth;
        etai = viewray.eta;
        etat = (*mtiter).eta;
        //cout <<"reflect etai=" <<etai <<endl;
        //cout <<"reflect etat=" <<etat <<endl;
        F0 = pow((etat - etai) / (etat + etai), 2);
        Fr = F0 + (1.0f - F0) * pow(1.0f - N.dot(V), 5);
        RayType refleray;
        RayPayType refleraypay;
        vec3d R;
        R = (2.0f * N.dot(V) * N - V).normalize();
        refleray.rdir = R;
        refleray.rorg = intersect + bias * refleray.rdir;
        refleray.eta = viewray.eta;
        Trace_Ray(refleray, scene, sphere, triangle, refleraypay);
        //cout <<"refle inter object idx=" <<refleraypay.object_idx <<endl;
        //cout <<"refle inter object type=" <<refleraypay.object_type <<endl;
        //cout <<"refle inter object subidx=" <<refleraypay.subobject_idx <<endl;
        if(refleraypay.object_idx > -1)
        {
            refledepth++;
            //cout <<"reflect depth=" <<refledepth <<endl;
            reflecolor = Shade_Ray(refleray, refleraypay, scene, sphere, triangle, material, texture, light, img, refledepth);
        }
        I += Fr * reflecolor;
        
        vec3d refracolor(0.0f);
        int refradepth;
        refradepth = depth;
        RayType refraray;
        RayPayType refraraypay;
        vec3d T;
        float alpha = 0.0f, t = 0.0f;
        alpha = (*mtiter).alpha;
        t = raypay.distance;
        if(N.dot(V) > 0.0f)
        {
            T = (N * (-1) * sqrt(1.0f - pow(etai / etat, 2) * (1.0f - pow(N.dot(V), 2))) + etai / etat * (N.dot(V) * N - I)).normalize();
            refraray.rdir = T;
            refraray.rorg = intersect + bias * refraray.rdir;
            refraray.eta = etat;
            Trace_Ray(refraray, scene, sphere, triangle, refraraypay);
            if(refraraypay.object_idx > -1)
            {
                //cout <<"refract inter object idx=" <<refraraypay.object_idx <<endl;
                //cout <<"refract inter object type=" <<refraraypay.object_type <<endl;
                //cout <<"refract inter object subidx=" <<refraraypay.subobject_idx <<endl;
                refradepth++;
                //cout <<"refract in depth=" <<refradepth <<endl;
                refracolor = Shade_Ray(refraray, refraraypay, scene, sphere, triangle, material, texture, light, img, refradepth);
            }
            F0 = pow((etat - etai) / (etat + etai), 2);
            Fr = F0 + (1.0f - F0) * pow(1.0f - N.dot(V), 5);
        }
        else if(N.dot(V) < 0.0f)
        {
            //cout <<"cosi=" <<N.dot(V) <<endl;
            //cout <<"refract out etai=" <<etai <<endl;
            //cout <<"refract out etat=" <<etat <<endl;
            etat = 1.0f;
            N = N * (-1);
            float k;
            k = 1.0f - pow(etai / etat, 2) * (1.0f - pow(N.dot(V), 2));
            if(k > 0.0f || k == 0.0f)
            {
                //cout <<"refract out" <<endl;
                T = (N * (-1) * sqrt(1.0f - pow(etat / etai, 2) * (1.0f - pow(N.dot(V), 2))) + (etat / etai) * (N.dot(V) * N - I)).normalize();
                refraray.rdir = T;
                refraray.eta = etat;
                Trace_Ray(refraray, scene, sphere, triangle, refraraypay);
                if(refraraypay.object_idx > -1)
                {
                    refradepth++;
                    //cout <<"refract out depth=" <<refradepth <<endl;
                    refracolor = Shade_Ray(refraray, refraraypay, scene, sphere, triangle, material, texture, light, img, refradepth);
                }
                F0 = pow((etat - etai) / (etat + etai), 2);
                Fr = F0 + (1.0f - F0) * pow(1.0f - N.dot(V), 5);
            }
            else
            {
                Fr = 1.0f;
                //cout <<"total internel reflection" <<endl;
            }
        }
        else
        {
            Fr = 1.0f;
            //cout <<"tangent reflection" <<endl;
        }
        //I += (1.0f - Fr) * exp((-1.0f) * alpha * t) * refracolor;
        I += (1.0f - Fr) * (1.0f - alpha) * refracolor;
    }
    
    if(I.x > 1.0f)
    {
        I.x = 1.0f;
    }
    if(I.y > 1.0f)
    {
        I.y = 1.0f;
    }
    if(I.z > 1.0f)
    {
        I.z = 1.0f;
    }
    img = I;
    return img;
}

void Sphere_Texture(const SphereType sphere, const vec3d intersect, float & u, float & v)
{
    float Nx, Ny, Nz;
    float phi, theta;
    float pi = 4.0f * atan(1.0f);
    Nx = (intersect.x - sphere.sorg.x) / sphere.r;
    Ny = (intersect.y - sphere.sorg.y) / sphere.r;
    Nz = (intersect.z - sphere.sorg.z) / sphere.r;
    phi = acos(-Ny);
    theta = atan2(-Nx, -Nz);
    v = phi / pi;
    if(theta < 0.0f)
    {
        theta += 2.0f * pi;
    }
    u = theta / (2.0f * pi);
}

void Triangle_Texture(const TriangleType triangle, float & u, float & v)
{
    float u0, v0, u1, v1, u2, v2;
    u0 = triangle.vertexture1.x;
    v0 = triangle.vertexture1.y;
    u1 = triangle.vertexture2.x;
    v1 = triangle.vertexture2.y;
    u2 = triangle.vertexture3.x;
    v2 = triangle.vertexture3.y;
    
    float alpha, beta, gamma;
    alpha = triangle.alpha;
    beta = triangle.beta;
    gamma = triangle.gamma;
    u = alpha * u0 + beta * u1 + gamma * u2;
    v = alpha * v0 + beta * v1 + gamma * v2;
}
