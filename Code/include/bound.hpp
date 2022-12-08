#ifndef BOUND_H
#define BOUND_H
#include <vecmath.h>
#include <vector>
#include "ray.hpp"
using std::vector;

class AABB
{
public:
    AABB();
    AABB(Vector3f &low, Vector3f &high);
    bool intersect(const Ray &r, float &collision_t_min, float t_min);
    void update(const Vector3f &vec);
    void init()
    {
        min_b = Vector3f::INFVEC;
        max_b = (-1.0) * Vector3f::INFVEC;
    }
    Vector3f get_min() const { return min_b; }
    Vector3f get_max() const { return max_b; }
    void set(const Vector3f &lo, const Vector3f &hi)
    {
        min_b = lo;
        max_b = hi;
    }
    Vector3f *get_8_vertices() //获取八个顶点位置
    {
        Vector3f *verts = new Vector3f[8];
        for (int i = 0; i < 8; i++)
        {
            if (i % 2 == 0)
                verts[i][2] = min_b[2];
            else
                verts[i][2] = max_b[2];
            if (i % 4 < 2)
                verts[i][1] = min_b[1];
            else
                verts[i][1] = max_b[1];
            if (i < 4)
                verts[i][0] = min_b[0];
            else
                verts[i][0] = max_b[0];
        }
        return verts;
    }

    bool no_point_intsc(const Vector3f &p, float sensitive_radius) //光是否能给某个hitpoint子树送光的快速判定
    {
        float dist = 0;
        for (int i = 0; i < 3; i++)
        {
            if (p[i] > max_b[i])
                dist += (p[i] - max_b[i]) * (p[i] - max_b[i]);
            else if (p[i] < min_b[i])
                dist += (p[i] - min_b[i]) * (p[i] - min_b[i]);
        }
        if (dist > sensitive_radius)
            return true;
        return false;
    }

    void split(Vector3f &left_max, Vector3f &right_min, int dim) // object seg tree建树中，子树参数的求解
    {
        left_max = max_b;
        right_min = min_b;
        left_max[dim] = (min_b[dim] + max_b[dim]) / 2.0;
        right_min[dim] = (min_b[dim] + max_b[dim]) / 2.0;
    }

    bool involved(Vector3f &lo, Vector3f &hi) //判定某一个包围盒是否被这个包围盒包围
    {
        for (int i = 0; i < 3; i++)
        {
            bool tmp1 = (lo[i] < max_b[i]);
            bool tmp2 = (lo[i] == hi[i] && lo[i] == max_b[i]);
            if (!(tmp1 || tmp2))
                return false;
            tmp1 = (hi[i] > min_b[i]);
            tmp2 = (lo[i] == hi[i] && hi[i] == min_b[i]);
            if (!(tmp1 || tmp2))
                return false;
        }
        return true;
    }

private:
    Vector3f min_b;
    Vector3f max_b;
    void refine(float &curr_min, float &curr_max);
    inline float perturbe(float x) //防止出现除以零的情况
    {
        if (x == 0.0)
            x = x + 1e-10;
        return x;
    }
};

void AABB::update(const Vector3f &vec) //更新aabb bounding box参数
{
    for (int i = 0; i < 3; i++)
    {
        min_b[i] = min_b[i] < vec[i] ? min_b[i] : vec[i];
        max_b[i] = max_b[i] > vec[i] ? max_b[i] : vec[i];
    }
}

void AABB::refine(float &curr_min, float &curr_max)
{
    if (curr_min > curr_max) //如果顺序反了就交换一下
    {
        float tmp = curr_max;
        curr_max = curr_min;
        curr_min = tmp;
    }
    return;
}
AABB::AABB(Vector3f &low, Vector3f &high)
{
    min_b = low;
    max_b = high;
}

AABB::AABB()
{
    min_b = Vector3f::INFVEC;
    max_b = (-1.0) * Vector3f::INFVEC;
}

bool AABB::intersect(const Ray &r, float &collision_t_min, float t_min) //求解视线是否与bounding box相交
{
    // collision_t_min是一个记录相交位置最小值的功能，可以用于迅速排除一些不可能相交的位置
    float t_x_min, t_x_max, t_y_min, t_y_max, t_z_min, t_z_max;
    float collision_min, collision_max;

    //(x-ox)/dx
    t_x_min = (min_b.x() - r.getOrigin().x()) / perturbe(r.getDirection().x());
    t_x_max = (max_b.x() - r.getOrigin().x()) / perturbe(r.getDirection().x());
    refine(t_x_min, t_x_max);
    collision_min = t_x_min;
    collision_max = t_x_max;

    t_y_min = (min_b.y() - r.getOrigin().y()) / perturbe(r.getDirection().y());
    t_y_max = (max_b.y() - r.getOrigin().y()) / perturbe(r.getDirection().y());
    refine(t_y_min, t_y_max);

    if (t_x_min > t_y_max || t_y_min > t_x_max)
        return false;

    collision_min = collision_min < t_y_min ? t_y_min : collision_min;
    collision_max = collision_max > t_y_max ? t_y_max : collision_max;

    t_z_min = (min_b.z() - r.getOrigin().z()) / perturbe(r.getDirection().z());
    t_z_max = (max_b.z() - r.getOrigin().z()) / perturbe(r.getDirection().z());
    refine(t_z_min, t_z_max);

    if (collision_min > t_z_max || t_z_min > collision_max)
        return false;

    collision_t_min = collision_min < t_z_min ? t_z_min : collision_min;
    float collision_t_max = collision_max > t_z_max ? t_z_max : collision_max;

    if (collision_t_max < t_min) //这样几乎不会出现负向相交的盒子，不知道这样加会不会出错
        return false;
    return true;
}

#endif // !BOUND_H