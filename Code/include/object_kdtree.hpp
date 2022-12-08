#ifndef OBJECTKDTREE_H
#define OBJECTKDTREE_H
#include <vecmath.h>

/*参考开源代码
https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project */

#include <map>
#include <vector>

#include "bound.hpp"
#include "hit.hpp"
#include "object3d.hpp"

#define MAX_DEPTH 24
#define MAX_FACES 128

#ifndef INF
#define INF 0x7fffffff;
#endif

using std::map;
using std::vector;
class ObjectSegTreeNode
{
public:
    AABB aabb;
    vector<Object3D *> *faces;
    ObjectSegTreeNode *lc, *rc;
    int l, r;
    bool involved(Object3D *face)
    {
        Vector3f faceMin = face->min();
        Vector3f faceMax = face->max();
        return aabb.involved(faceMin, faceMax);
    }
};

class ObjectSegTree
{
    int n;
    Vector3f **vertices;
    ObjectSegTreeNode *create_node(vector<Object3D *> *faces, int depth, int d,
                                   const Vector3f &min, const Vector3f &max)
    {
        ObjectSegTreeNode *p = new ObjectSegTreeNode;
        p->aabb.set(min, max);
        Vector3f maxL, minR;
        p->aabb.split(maxL, minR, d);
        p->faces = new vector<Object3D *>;
        for (int i = 0; i < (*faces).size(); i++)
            if (p->involved((*faces)[i]))
                p->faces->push_back((*faces)[i]);

        if (p->faces->size() > MAX_FACES && depth < MAX_DEPTH) //深度或者面数突破最大限制
        {
            p->lc = create_node(p->faces, depth + 1, (d + 1) % 3, min, maxL);
            p->rc = create_node(p->faces, depth + 1, (d + 1) % 3, minR, max);

            vector<Object3D *> *faceL = p->lc->faces, *faceR = p->rc->faces;
            map<Object3D *, int> cnt;
            for (auto face : *faceL)
                cnt[face]++;
            for (auto face : *faceR)
                cnt[face]++;
            p->lc->faces = new vector<Object3D *>;
            p->rc->faces = new vector<Object3D *>;
            p->faces->clear();
            for (auto face : *faceL)
                if (cnt[face] == 1)
                    p->lc->faces->push_back(face);
                else
                    p->faces->push_back(face);
            for (auto face : *faceR)
                if (cnt[face] == 1)
                    p->rc->faces->push_back(face);
        }
        else
            p->lc = p->rc = nullptr;
        return p;
    }

    void getFaces(ObjectSegTreeNode *p, vector<Object3D *> *faces)
    {
        p->l = faces->size();
        for (auto face : *(p->faces))
            faces->push_back(face);
        p->r = faces->size();
        if (p->lc)
            getFaces(p->lc, faces);
        if (p->rc)
            getFaces(p->rc, faces);
    }

public:
    ObjectSegTreeNode *root;
    vector<Object3D *> *faces;
    ObjectSegTree(vector<Object3D *> *faces)
    {
        Vector3f min = Vector3f(INF, INF, INF);
        Vector3f max = -min;
        for (auto face : *faces)
        {
            min = minE(min, face->min());
            max = maxE(max, face->max());
        }
        root = create_node(faces, 1, 0, min, max);
        this->faces = new vector<Object3D *>;
        getFaces(root, this->faces);
    }

    float cuboidIntersect(ObjectSegTreeNode *p, const Ray &ray) const
    {
        if (!p)
            return INF;
        float t = INF;
        p->aabb.intersect(ray, t, 0);
        return t;
    }

    bool intersect(const Ray &ray, Hit &hit) const
    {
        Object3D *nextFace = nullptr;
        return intersect(root, ray, nextFace, hit);
    }

    bool intersect(ObjectSegTreeNode *p, const Ray &ray, Object3D *&nextFace,
                   Hit &hit) const
    {
        bool flag = false;
        for (int i = 0; i < p->faces->size(); ++i)
            if ((*p->faces)[i]->intersect(ray, hit))
            {
                nextFace = (*p->faces)[i];
                flag = true;
            }
        float tl = cuboidIntersect(p->lc, ray),
              tr = cuboidIntersect(p->rc, ray);
        if (tl < tr) //逻辑冗余
        {
            if (hit.t <= tl)
                return flag;
            if (p->lc)
                flag |= intersect(p->lc, ray, nextFace, hit);
            if (hit.t <= tr)
                return flag;
            if (p->rc)
                flag |= intersect(p->rc, ray, nextFace, hit);
        }
        else
        {
            if (hit.t <= tr)
                return flag;
            if (p->rc)
                flag |= intersect(p->rc, ray, nextFace, hit);
            if (hit.t <= tl)
                return flag;
            if (p->lc)
                flag |= intersect(p->lc, ray, nextFace, hit);
        }
        return flag;
    }
};

#endif // !OBJECTKDTREE_H