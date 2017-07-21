//
//  main.cpp
//  ascii_image
//
//  Created by Alex on 9/11/16.
//  Copyright © 2016 Alex. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "vec3d.h"
#include "raytrace.h"

using namespace std;

int main(int argc, const char * argv[])
{
    string filename = "";
    cout <<"input file name:" <<endl;
    cin >>filename;
    
    char inputname[100] = "";  // store path of output file
    //sprintf(inputname, "/Users/alex/Desktop/%s.txt", filename.c_str());
    sprintf(inputname, "./%s.txt", filename.c_str());
    cout <<"input file:" <<inputname <<endl;
    ifstream fin(inputname);
    
    char outputname[100] = "";  // store path of output file
    //sprintf(outputname, "/Users/alex/Desktop/%s.ppm", filename.c_str());
    sprintf(outputname, "./%s.ppm", filename.c_str());
    ofstream fout(outputname);
    
    RayType ray;
    vec3d up;
    float fovv = -1.0f;
    bool parallel = false;
    int width = -1, height = -1;
    vec3d bcolor;
    vector<MaterialType> material;
    vector<SphereType> sphere;
    vector<LightType> light;
    vector<vec3d> vertex;
    vector<vec3d> vertexure;
    vector<vec3d> vertexnormal;
    vector<TriangleType> triangle;
    vector<string> texname;
    vector<TextureType> texture;
    int material_num = 0, sphere_num = 0, light_num = 0, vertex_num = 0, vertexure_num = 0, vertexnormal_num = 0, triangle_num = 0, texture_num = 0;
    int i, j, k;
    
    if(fin.is_open())
    {
        string s = "";
        while(getline(fin, s))  // read input file by line
        {
            string s2(s);
            string s3(s);
            string s4(s);
            /*cout <<"s=" <<s <<endl;
             cout <<"s2=" <<s2 <<endl;
             cout <<"s3=" <<s3 <<endl;
             cout <<"s4=" <<s4 <<endl;*/
            stringstream line(s);
            stringstream line2(s2);
            stringstream line3(s3);
            stringstream line4(s4);
            string var = "";
            line >>var;
            
            // set parameter from the input file to variable
            if(var == "eye")
            {
                line >>ray.rorg.x >>ray.rorg.y >>ray.rorg.z;
                /*cout <<"eyex=" <<ray.rorg.x <<endl;
                 cout <<"eyey=" <<ray.rorg.y <<endl;
                 cout <<"eyez=" <<ray.rorg.z <<endl;*/
            }
            else if(var == "viewdir")
            {
                line >>ray.rdir.x >>ray.rdir.y >>ray.rdir.z;
                /*cout <<"viewx=" <<ray.rdir.x <<endl;
                 cout <<"viewy=" <<ray.rdir.y <<endl;
                 cout <<"viewz=" <<ray.rdir.z <<endl;*/
            }
            else if(var == "updir")
            {
                line >>up.x >>up.y >>up.z;
                /*cout <<"upx=" <<up.x <<endl;
                 cout <<"upy=" <<up.y <<endl;
                 cout <<"upz=" <<up.z <<endl;*/
            }
            else if(var == "fovv")
            {
                line >>fovv;
                //cout <<"fovv=" <<fovv <<endl;
            }
            else if(var =="parallel")
            {
                parallel = true;
            }
            else if(var == "imsize")
            {
                line >>width >>height;
                //cout <<"width="<<width <<endl;
                //cout <<"height=" <<height <<endl;
            }
            else if(var == "bkgcolor")
            {
                line >>bcolor.x >>bcolor.y >>bcolor.z;
                /*cout <<"bcolorx=" <<bcolor.x <<endl;
                 cout <<"bcolory=" <<bcolor.y <<endl;
                 cout <<"bcolorz=" <<bcolor.z <<endl;*/
            }
            else if(var == "mtlcolor")
            {
                MaterialType mtl1;
                mtl1.idx = material_num;
                line >>mtl1.dcolor.x >>mtl1.dcolor.y >>mtl1.dcolor.z >>mtl1.scolor.x >>mtl1.scolor.y >>mtl1.scolor.z >>mtl1.ka >>mtl1.kd >>mtl1.ks >>mtl1.n >>mtl1.alpha >>mtl1.eta;
                material.push_back(mtl1);
                /*cout <<"dcolorx=" <<mtl1.dcolor.x <<endl;
                 cout <<"dcolory=" <<mtl1.dcolor.y <<endl;
                 cout <<"dcolorz=" <<mtl1.dcolor.z <<endl;
                 cout <<"scolorx=" <<mtl1.scolor.x <<endl;
                 cout <<"scolory=" <<mtl1.scolor.y <<endl;
                 cout <<"scolorz=" <<mtl1.scolor.z <<endl;
                 cout <<"ka=" <<mtl1.ka <<endl;
                 cout <<"kd=" <<mtl1.kd <<endl;
                 cout <<"ks=" <<mtl1.ks <<endl;
                 cout <<"n=" <<mtl1.n <<endl;*/
                material_num++;
            }
            else if(var == "texture")
            {
                string texname1;
                line >>texname1;
                texname.push_back(texname1);
                //cout <<"texture name=" <<texname1 <<endl;
                texture_num++;
                //cout <<"texture_num=" <<texture_num <<endl;
            }
            else if(var == "sphere")
            {
                SphereType sph1;
                sph1.idx = sphere_num;
                line >>sph1.sorg.x >>sph1.sorg.y >>sph1.sorg.z >>sph1.r;
                /*cout <<"sphx=" <<sph1.sorg.x <<endl;
                 cout <<"sphy=" <<sph1.sorg.y <<endl;
                 cout <<"sphz=" <<sph1.sorg.z <<endl;
                 cout <<"sphr=" <<sph1.r <<endl;*/
                sph1.material = material_num - 1;
                if(texture_num > 0)
                {
                    sph1.texture = texture_num - 1;
                }
                else
                {
                    sph1.texture = -1;
                }
                sphere.push_back(sph1);
                cout <<"sphere material=" <<sph1.material <<endl;
                cout <<"sphere texture=" <<sph1.texture <<endl;
                sphere_num++;
                //cout <<"sphere num=" <<sphere_num <<endl;
            }
            else if(var == "light")
            {
                LightType lgt1;
                lgt1.idx = light_num;
                line >>lgt1.lorg.x >>lgt1.lorg.y >>lgt1.lorg.z >>lgt1.w >>lgt1.lcolor.x >>lgt1.lcolor.y >>lgt1.lcolor.z;
                light.push_back(lgt1);
                /*cout <<"lorgx=" <<lgt1.lorg.x <<endl;
                 cout <<"lorgy=" <<lgt1.lorg.y <<endl;
                 cout <<"lorgz=" <<lgt1.lorg.z <<endl;
                 cout <<"lw=" <<lgt1.w <<endl;
                 cout <<"lcolorx=" <<lgt1.lcolor.x <<endl;
                 cout <<"lcolory=" <<lgt1.lcolor.y <<endl;
                 cout <<"lcolorz=" <<lgt1.lcolor.z <<endl;*/
                light_num++;
            }
            else if(var == "v")
            {
                vec3d vertex1;
                line >>vertex1.x >>vertex1.y >>vertex1.z;
                vertex.push_back(vertex1);
                /*cout <<"v.x=" <<vertex1.x <<endl;
                 cout <<"v.y=" <<vertex1.y <<endl;
                 cout <<"v.z=" <<vertex1.z <<endl;*/
                vertex_num++;
                //cout <<"vertex num=" <<vertex_num <<endl;
            }
            else if(var == "vt")
            {
                vec3d vertexure1;
                line >>vertexure1.x >>vertexure1.y;
                vertexure1.z = 0.0f;
                vertexure.push_back(vertexure1);
                //cout <<"vt.x=" <<vertexure1.x <<endl;
                //cout <<"vt.y=" <<vertexure1.y <<endl;
                vertexure_num++;
            }
            else if(var == "vn")
            {
                vec3d vertexnormal1;
                line >>vertexnormal1.x >>vertexnormal1.y >>vertexnormal1.z;
                vertexnormal.push_back(vertexnormal1);
                /*cout <<"vn.x=" <<vertexnormal1.x <<endl;
                 cout <<"vn.y=" <<vertexnormal1.y <<endl;
                 cout <<"vn.z=" <<vertexnormal1.z <<endl;*/
                vertexnormal_num++;
            }
            else if(var == "f")
            {
                TriangleType triangle1;
                triangle1.idx = triangle_num;
                char symbol;
                line >>triangle1.v1 >>symbol;
                //cout <<"symbol=" <<symbol <<endl;
                if(symbol == '/')
                {
                    line >>symbol;
                    //cout <<"symbol2=" <<symbol <<endl;
                    if(symbol == '/')
                    {
                        line >>triangle1.vn1 >>triangle1.v2 >>symbol >>symbol >>triangle1.vn2 >>triangle1.v3 >>symbol >>symbol >>triangle1.vn3;
                        triangle1.texture = -1;
                        triangle1.smooth = true;
                        /*cout <<"v1=" <<triangle1.v1 <<endl;
                         cout <<"v2=" <<triangle1.v2 <<endl;
                         cout <<"v3=" <<triangle1.v3 <<endl;
                         cout <<"vn1=" <<triangle1.vn1 <<endl;
                         cout <<"vn2=" <<triangle1.vn2 <<endl;
                         cout <<"vn3=" <<triangle1.vn3 <<endl;*/
                    }
                    else
                    {
                        line2 >>var >>triangle1.v1 >>symbol >>triangle1.vt1 >>symbol;
                        //cout <<"symbol3=" <<symbol <<endl;
                        if(symbol == '/')
                        {
                            line2 >>triangle1.vn1 >>triangle1.v2 >>symbol >>triangle1.vt2 >>symbol >>triangle1.vn2 >>triangle1.v3 >>symbol >>triangle1.vt3 >>symbol >>triangle1.vn3;
                            triangle1.texture = texture_num - 1;
                            triangle1.smooth = true;
                            /*cout <<"v1=" <<triangle1.v1 <<endl;
                             cout <<"v2=" <<triangle1.v2 <<endl;
                             cout <<"v3=" <<triangle1.v3 <<endl;
                             cout <<"vt1=" <<triangle1.vt1 <<endl;
                             cout <<"vt2=" <<triangle1.vt2 <<endl;
                             cout <<"vt3=" <<triangle1.vt3 <<endl;
                             cout <<"vn1=" <<triangle1.vn1 <<endl;
                             cout <<"vn2=" <<triangle1.vn2 <<endl;
                             cout <<"vn3=" <<triangle1.vn3 <<endl;*/
                        }
                        else
                        {
                            line3 >>var >>triangle1.v1 >>symbol >>triangle1.vt1 >>triangle1.v2 >>symbol >>triangle1.vt2 >>triangle1.v3 >>symbol >>triangle1.vt3;
                            triangle1.texture = texture_num - 1;
                            triangle1.smooth = false;
                            /*cout <<"v1=" <<triangle1.v1 <<endl;
                             cout <<"v2=" <<triangle1.v2 <<endl;
                             cout <<"v3=" <<triangle1.v3 <<endl;
                             cout <<"vt1=" <<triangle1.vt1 <<endl;
                             cout <<"vt2=" <<triangle1.vt2 <<endl;
                             cout <<"vt3=" <<triangle1.vt3 <<endl;*/
                        }
                    }
                }
                else
                {
                    line4 >>var >>triangle1.v1 >>triangle1.v2 >>triangle1.v3;
                    triangle1.texture = -1;
                    triangle1.smooth = false;
                    /*cout <<"v1=" <<triangle1.v1 <<endl;
                     cout <<"v2=" <<triangle1.v2 <<endl;
                     cout <<"v3=" <<triangle1.v3 <<endl;*/
                }
                triangle1.material = material_num - 1;
                triangle.push_back(triangle1);
                cout <<"triangle material=" <<triangle1.material <<endl;
                cout <<"triangle texture=" <<triangle1.texture <<endl;
                cout <<"triangle smooth=" <<triangle1.smooth <<endl;
                triangle_num++;
                //cout <<"triangle num=" <<triangle_num <<endl;
            }
            else if(var == "" || var == "#"){}
            else
            {
                cerr <<"Can not recognize the input!" <<endl;
                break;
            }
        }
    }
    // print error when the input file can not be found
    else
    {
        cerr <<"Can not find file " <<filename <<" !" <<endl;
    }
    fin.close();
    
    // when the input value is valid, generate image to output file
    if(fovv < 180.0f && width > 0 && height > 0)
    {
        vector< vector<vec3d> > img(1000, vector<vec3d>(1000, bcolor));
        
        if(triangle_num > 0)
        {
            vector<vec3d>::iterator viter = vertex.begin();
            vector<TriangleType>::iterator triter;
            for(triter = triangle.begin(); triter != triangle.end(); triter++)
            {
                (*triter).vertex1 = *(viter + (*triter).v1 - 1);
                (*triter).vertex2 = *(viter + (*triter).v2 - 1);
                (*triter).vertex3 = *(viter + (*triter).v3 - 1);
                if(vertexure_num > 0)
                {
                    vector<vec3d>::iterator vtxureiter = vertexure.begin();
                    (*triter).vertexture1 = *(vtxureiter + (*triter).vt1 - 1);
                    (*triter).vertexture2 = *(vtxureiter + (*triter).vt2 - 1);
                    (*triter).vertexture3 = *(vtxureiter + (*triter).vt3 - 1);
                }
                if(vertexnormal_num > 0)
                {
                    vector<vec3d>::iterator vtxnormaliter = vertexnormal.begin();
                    (*triter).normal1 = *(vtxnormaliter + (*triter).vn1 - 1);
                    (*triter).normal2 = *(vtxnormaliter + (*triter).vn2 - 1);
                    (*triter).normal3 = *(vtxnormaliter + (*triter).vn3 - 1);
                }
                /*cout <<"triangle v1=" <<(*triter).vertex1 <<endl;
                 cout <<"triangle v2=" <<(*triter).vertex2 <<endl;
                 cout <<"triangle v3=" <<(*triter).vertex3 <<endl;
                 cout <<"triangle vt1=" <<(*triter).vertexture1 <<endl;
                 cout <<"triangle vt2=" <<(*triter).vertexture2 <<endl;
                 cout <<"triangle vt3=" <<(*triter).vertexture3 <<endl;
                 cout <<"triangle vn1=" <<(*triter).normal1 <<endl;
                 cout <<"triangle vn2=" <<(*triter).normal2 <<endl;
                 cout <<"triangle vn3=" <<(*triter).normal3 <<endl;*/
            }
        }
        
        if(texture_num > 0)
        {
            vector<string>::iterator texnaiter;
            for(texnaiter = texname.begin(); texnaiter != texname.end(); texnaiter++)
            {
                TextureType texture1;
                char texturename[100] = "";  // store path of output file
                //sprintf(texturename, "/Users/alex/Desktop/%s", (*texnaiter).c_str());
                sprintf(texturename, "./%s", (*texnaiter).c_str());
                cout <<"texture file:" <<texturename <<endl;
                ifstream texturein(texturename);
                vec3d texturecolor[10000];
                if(texturein.is_open())
                {
                    string s = "";
                    i = 0;
                    while(getline(texturein, s))  // read input file by line
                    {
                        stringstream line(s);
                        string var = "";
                        line >>var;
                        if(var == "P3")
                        {
                            line >>texture1.width >>texture1.height;
                        }
                        else
                        {
                            float x, y, z;
                            stringstream ss(var);
                            ss >>x;
                            line >>y >>z;
                            texturecolor[i] = vec3d(x, y, z);
                            i++;
                        }
                    }
                }
                else
                {
                    cerr <<"Can not find texture file " <<*texnaiter <<" !" <<endl;
                }
                texturein.close();
                
                k = 0;
                for(i = 0; i < texture1.height; i++)
                {
                    for(j = 0; j < texture1.width; j++)
                    {
                        texture1.texcolor[i][j] = texturecolor[k];
                        k++;
                    }
                }
                texture.push_back(texture1);
            }
        }
        
        vector<SceneType> scene;
        int scene_num;
        i = 0;
        if(triangle_num > 0)
        {
            vector<TriangleType>::iterator triter;
            for(triter = triangle.begin(); triter != triangle.end(); triter++)
            {
                SceneType scene1;
                scene1.object_idx = i;
                scene1.object_type = 1;
                scene1.subobject_idx = (*triter).idx;
                scene1.material_idx = (*triter).material;
                scene1.texture_idx = (*triter).texture;
                scene.push_back(scene1);
                i++;
            }
        }
        if(sphere_num > 0)
        {
            vector<SphereType>::iterator spiter;
            for(spiter = sphere.begin(); spiter != sphere.end(); spiter++)
            {
                SceneType scene1;
                scene1.object_idx = i;
                scene1.object_type = 0;
                scene1.subobject_idx = (*spiter).idx;
                scene1.material_idx = (*spiter).material;
                scene1.texture_idx = (*spiter).texture;
                scene.push_back(scene1);
                i++;
            }
        }
        scene_num = i;
        cout <<"scene num=" <<scene_num <<endl;
        
        for(i = 0; i < height; i++)
        {
            for(j = 0; j < width; j++)
            {
                RayType viewray;
                viewray.eta = 1.0f;
                RayPayType raypay;
                Get_Ray(ray, up, i, j, width, height, fovv, parallel, viewray);
                Trace_Ray(viewray, scene, sphere, triangle, raypay);
                if(raypay.object_idx > -1)
                {
                    //cout <<"inter idx=" <<raypay.object_idx <<endl;
                    //cout <<"inter type=" <<raypay.object_type <<endl;
                    //cout <<"inter subidx=" <<raypay.subobject_idx <<endl;
                    int depth = 1;
                    img[i][j] = Shade_Ray(viewray, raypay, scene, sphere, triangle, material, texture, light, img[i][j], depth);
                }
            }
        }
        
        fout <<"P3" <<endl;
        fout <<"# ascii image" <<endl;
        fout <<width <<"\t" <<height <<endl;
        fout <<255 <<endl;
        for(i = 0; i< height; i++)
        {
            for(j = 0; j < width; j++)
            {
                img[i][j] *= 255.0f;
                fout <<int(img[i][j].x) <<"\t" <<int(img[i][j].y) <<"\t" <<int(img[i][j].z) <<"\t" <<endl;
            }
        }
    }
    // print error when the input value is invalid
    else
    {
        cerr <<"Input value is invalid, can not generate image!" <<endl;
    }
    
    fout.close();
    return 0;
}
