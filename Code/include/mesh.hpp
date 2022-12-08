#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <vecmath.h>

#include "object3d.hpp"
#include "object_kdtree.hpp"

#include "bound.hpp"
#include "ray.hpp"
#include "triangle.hpp"
#include "utils.hpp"

/* 参考 https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project */

class Mesh : public Object3D
{

public:
    Mesh(const char *filename, Material *m) : Object3D(m)
    {
        set_normal = false;
        set_uv = false;
        std::ifstream f;
        f.open(filename);
        if (!f.is_open())
        {
            std::cout << "Cannot open " << filename << "\n";
            return;
        }
        std::string line;
        std::string vTok("v");
        std::string fTok("f");
        std::string normalTok("vn"); // normal信息
        std::string texTok("vt");    //纹理信息
        std::string bslash("/"), space(" ");
        std::string tok;
        int texID;

        while (true)
        {
            std::getline(f, line);
            if (f.eof())
            {
                break;
            }
            if (line.size() < 3)
            {
                continue;
            }
            if (line.at(0) == '#')
            {
                continue;
            }
            std::stringstream ss(line);
            ss >> tok;

            if (tok == vTok)
            {
                Vector3f vec;
                ss >> vec[0] >> vec[1] >> vec[2];
                aabb.update(vec); //
                v.push_back(vec);
            }
            else if (tok == fTok)
            {
                TriangleIndex vId, tId, nId;
                for (int i = 0; i < 3; ++i)
                {
                    std::string str;
                    ss >> str;
                    std::vector<std::string> id = split(str, bslash);
                    vId[i] = atoi(id[0].c_str()) - 1;
                    if (id.size() > 1)
                    {
                        tId[i] = atoi(id[1].c_str()) - 1;
                    }
                    if (id.size() > 2)
                    {
                        nId[i] = atoi(id[2].c_str()) - 1;
                    }
                }
                vIdx.push_back(vId);
                tIdx.push_back(tId);
                nIdx.push_back(nId);
            }
            else if (tok == texTok)
            {
                Vector2f texcoord;
                ss >> texcoord[0];
                ss >> texcoord[1];
                vt.push_back(texcoord);
            }
            else if (tok == normalTok)
            {
                Vector3f vec;
                ss >> vec[0] >> vec[1] >> vec[2];
                vn.push_back(vec);
            }
        }
        f.close();
        match();

        kdTree = new ObjectSegTree(&triangles);
        vIdx.clear();
        tIdx.clear();
        nIdx.clear();
        v.clear();
        vn.clear();
        vt.clear();
    }

    ~Mesh()
    {
        for (int i = 0; i < triangles.size(); ++i)
            delete triangles[i];
        delete kdTree;
    }

    struct TriangleIndex
    {
        TriangleIndex()
        {
            x[0] = -1;
            x[1] = -1;
            x[2] = -1;
        }
        int &operator[](const int i) { return x[i]; }
        int x[3]{};
    };

    vector<Object3D *> triangles;
    bool intersect(const Ray &r, Hit &h) override
    {
        float tb;
        if (!aabb.intersect(r, tb, 0))
            return false;
        if (tb > h.getT())
            return false;
        return kdTree->intersect(r, h);
    }

    Vector3f min() const override { return aabb.get_min(); }
    Vector3f max() const override { return aabb.get_max(); }
    Ray randomRay(int axis = -1, long long int seed = 0) const override
    {
        return Ray(Vector3f::ZERO, Vector3f::ZERO);
    }

private:
    AABB aabb;
    ObjectSegTree *kdTree;
    bool set_normal, set_uv;

    std::vector<TriangleIndex> vIdx, tIdx, nIdx;
    std::vector<Vector3f> v, vn;
    std::vector<Vector2f> vt;

    void match()
    {
        for (int triId = 0; triId < vIdx.size(); triId++)
        {
            TriangleIndex &v_id = vIdx[triId];
            TriangleIndex &uv_id = tIdx[triId];
            TriangleIndex &normal_id = nIdx[triId];
            triangles.push_back((Object3D *)new Triangle(
                v[v_id[0]], v[v_id[1]], v[v_id[2]], material));
            if (tIdx.size() != 0 && uv_id[0] != -1)
                ((Triangle *)triangles.back())->set_uvs(vt[uv_id[0]], vt[uv_id[1]], vt[uv_id[2]]);
            if (nIdx.size() != 0 && normal_id[0] != -1)
                ((Triangle *)triangles.back())->set_normals(vn[normal_id[0]], vn[normal_id[1]], vn[normal_id[2]]);
        }
    }
};

#endif
