#ifndef HIT_KDTREE_H
#define HIT_KDTREE_H

#include "bound.hpp"
#include "hit.hpp"

#ifndef ING
#define INF 0x7fffffff;
#endif

#include <algorithm>
using namespace std;

#define X_DIVIDE 0
#define Y_DIVIDE 1
#define Z_DIVIDE 2
#define ALPHA 0.7

class HitKDTreeNode
{
public:
    Hit *hit;
    AABB aabb;
    float radius_max; //这是可以捕捉光子的最大距离范围
    HitKDTreeNode *lc, *rc;
    ~HitKDTreeNode()
    {
        delete lc;
        delete rc;
    }
};

//参考开源代码：https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project */

class HitKDTree
{
    int n;
    Hit **hits;
    HitKDTreeNode *create_node(int l, int r, int d)
    {
        HitKDTreeNode *p = new HitKDTreeNode;
        p->rc = nullptr;
        p->lc = nullptr;
        p->radius_max = -INF;
        for (int i = l; i <= r; ++i)
        {
            p->aabb.update(hits[i]->p);                                                        //更新界标
            p->radius_max = (p->radius_max) > (hits[i]->r2) ? (p->radius_max) : (hits[i]->r2); //更新敏感范围
        }

        int mid = l + r >> 1;
        sort_hits(l, r, d);
        p->hit = hits[mid];

        if (l < mid)
            p->lc = create_node(l, mid - 1, (d + 1) % 3); //左子树存在
        if (mid < r)
            p->rc = create_node(mid + 1, r, (d + 1) % 3); //右子树存在
        return p;
    }

public:
    HitKDTreeNode *root;
    HitKDTree(vector<Hit *> *hits)
    {
        n = hits->size();
        this->hits = new Hit *[n];
        for (int i = 0; i < n; ++i)
            this->hits[i] = (*hits)[i];
        root = create_node(0, n - 1, 0);
    }
    ~HitKDTree()
    {
        if (!root)
            return;
        delete (root);
        delete[] hits;
    }

    void sort_hits(int start, int end, int mode)
    {
        int mid = (start + end) >> 1;
        if (mode == X_DIVIDE)
            nth_element(hits + start, hits + mid, hits + end + 1, cmpHitX);
        else if (mode == Y_DIVIDE)
            nth_element(hits + start, hits + mid, hits + end + 1, cmpHitY);
        else if (mode == Z_DIVIDE)
            nth_element(hits + start, hits + mid, hits + end + 1, cmpHitZ);
        else
        {
            std::cout << "Invalid division type!" << endl;
            exit(-1);
        }
    }

    void update(HitKDTreeNode *p, const Vector3f &photon,
                const Vector3f &attenuation, const Vector3f &d)
    {
        if (p == nullptr)
            return;

        //判定相交的情况
        if (p->aabb.no_point_intsc(photon, p->radius_max))
            return;

        if ((photon - p->hit->p).squaredLength() <= p->hit->r2) //如果确实有交点
        {
            Hit *hp = p->hit;
            hp->n++; //更新相交次数

            float factor = (ALPHA + hp->n * ALPHA) / ((hp->n) * ALPHA + 1.00);
            hp->flux = (hp->flux + hp->attenuation * attenuation) * factor; //更新光通量
            hp->r2 *= factor;                                               //更新感光半径
            p->radius_max = p->hit->r2;
        }
        if (p->lc)
        {
            update(p->lc, photon, attenuation, d);
            p->radius_max = (p->radius_max) > (p->lc->hit->r2) ? (p->radius_max) : (p->lc->hit->r2);
        }
        if (p->rc)
        {
            update(p->rc, photon, attenuation, d);
            p->radius_max = (p->radius_max) > (p->rc->hit->r2) ? (p->radius_max) : (p->rc->hit->r2);
        }
    }

    static bool cmpHitX(Hit *a, Hit *b) { return a->p.x() < b->p.x(); }

    static bool cmpHitY(Hit *a, Hit *b) { return a->p.y() < b->p.y(); }

    static bool cmpHitZ(Hit *a, Hit *b) { return a->p.z() < b->p.z(); }
};
#endif