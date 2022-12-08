#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include "bound.hpp"
#include <vecmath.h>
#include <cmath>
#include <assert.h>
using namespace std;

//这里加个aabb的快速判定

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D
{
public:
    Sphere()
    {
        _center = Vector3f(0.0, 0.0, 0.0);
        _radius = 1.0;
        // unit ball at the center
    }
    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material)
    {
        _center = center;
        _radius = radius;
        aabb.update(_center + radius);
        aabb.update(_center - radius);
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h) override
    {
        Vector3f r_origin = r.getOrigin();
        Vector3f r_dir = r.getDirection();
        Vector3f l = _center - r_origin;
        float l_len = l.length();
        float tp = l.dot(l, r.getDirection());
        float d;
        if (l_len * l_len - tp * tp > 0)
            d = sqrt(l_len * l_len - tp * tp);
        else
            d = 0.0;
        if (d > _radius)
            return false;
        float t_prime;
        if (_radius * _radius - d * d > 0)
            t_prime = sqrt(_radius * _radius - d * d);
        else
            t_prime = 0.0;

        if (l_len <= _radius) // inside
        {
            float curr_t = tp + t_prime;
            Vector3f itsct_point = r_origin + curr_t * r_dir;
            if ((curr_t >= 0) && (curr_t <= h.getT()))
            {
                Vector3f OP(r_origin + r_dir * curr_t - _center);
                Vector3f normal = OP.normalized();
                float u = 0.5 + atan2(normal.x(), normal.z()) / (2 * M_PI),
                      v = 0.5 - asin(normal.y()) / M_PI;
                h.set(curr_t, material, get_normal(normal, OP, u, v),
                      material->get_color(u, v), r_origin + r_dir * curr_t);
                return true;
            }
            else
                return false;
        }
        else // ouside
        {
            if (tp < 0)
                return false;
            float curr_t = tp - t_prime;
            Vector3f itsct_point = r.getOrigin() + curr_t * r.getDirection();
            if ((curr_t >= 0) && (curr_t <= h.getT()))
            {
                Vector3f OP(r_origin + r_dir * curr_t - _center);
                Vector3f normal = OP.normalized();
                float u = 0.5 + atan2(normal.x(), normal.z()) / (2 * M_PI),
                      v = 0.5 - asin(normal.y()) / M_PI;
                h.set(curr_t, material, get_normal(normal, OP, u, v),
                      material->get_color(u, v), r_origin + r_dir * curr_t);
                return true;
            }
            else
                return false;
        }
    }

    Vector3f get_normal(const Vector3f &n, const Vector3f &p, float u, float v)
    {
        Vector3f pu(-p.z(), 0, p.x());

        if (pu.squaredLength() < 1e-6)
            return n;

        float u_grad, v_grad, res;

        material->update_normal(u, v, u_grad, v_grad, res);

        if (fabs(res) < 1e-6)
            return n;

        Vector3f pv(cos(2 * M_PI * u) * p.y(), -sin(M_PI - M_PI * v) * _radius, sin(2 * M_PI * u) * p.y());
        Vector3f new_normal = Vector3f::cross(pu + n * u_grad / (2 * M_PI), pv + n * v_grad / M_PI);
        new_normal.normalize();
        return new_normal;
    }

    Ray randomRay(int axis = -1, long long int seed = 0) const override
    {
        float u = 2 * random(axis, seed) - 1, v = 2 * random(axis, seed) - 1;
        float r2 = u * u + v * v;
        while (r2 >= 1)
        {
            ++seed;
            u = 2 * random(axis, seed) - 1;
            v = 2 * random(axis, seed) - 1;
            r2 = u * u + v * v;
        }
        Vector3f dir(2 * u * sqrtf(1 - r2), 2 * v * sqrt(1 - r2), 1 - 2 * r2);
        dir.normalize();
        return Ray(_center + _radius * dir, dir);
    }

    Vector3f min() const override { return aabb.get_min(); }
    Vector3f max() const override { return aabb.get_max(); }

protected:
    Vector3f _center;
    float _radius;
    AABB aabb;
};

#endif
