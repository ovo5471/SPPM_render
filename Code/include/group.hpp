#pragma once

#include "hit.hpp"
#include "object3d.hpp"
#include "ray.hpp"
#include "object_kdtree.hpp"

#include <iostream>
#include <vector>

class Group
{
public:
    vector<Object3D *> illuminants;
    Group(const vector<Object3D *> &objs)
    {
        for (int i = 0; i < objs.size(); i++)
        {
            faces.push_back(objs[i]);
            if (objs[i]->material->emission != Vector3f::ZERO)
                illuminants.push_back(objs[i]);
        }
        kdTree = new ObjectSegTree(&faces); //创建最高一级的kd tree
    }
    ~Group() { delete kdTree; }

    // 对列表里所有物体都求一遍交点
    bool intersect(const Ray &r, Hit &h) { return kdTree->intersect(r, h); }

private:
    ObjectSegTree *kdTree;
    vector<Object3D *> faces;
};
