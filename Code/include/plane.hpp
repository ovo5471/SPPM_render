#ifndef PLANE_H
#define PLANE_H

#include <vecmath.h>

#include <cmath>

#include "object3d.hpp"
#ifndef INF
#define INF 0x7fffffff;
#endif

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D
{
public:
    Plane(const Vector3f &normal, float _d, Material *m) : Object3D(m)
    {
        _normal = normal.normalized();
        d = _d;
        ua = Vector3f::cross(Vector3f::UP, _normal);
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h) override
    {
        Vector3f o = r.getOrigin();
        Vector3f dir = r.getDirection();
        // dir.normalize();
        float cos = Vector3f::dot(_normal, dir);
        // 平行，这里原来有bug，原来的判定依据是cos大于0，但是大于0也是有可能相交的，导致有一些平面没有办法显示
        if (cos > -1e-6)
            return false;
        // if ((cos > -1e-6) && (cos < 1e-6))
        //     return false;
        float t = (d - Vector3f::dot(_normal, o)) / cos, u, v;

        if (t < 0 || t > h.getT())
            return false;

        Vector3f p(o + dir * t);
        v = p.y();
        v = v / 2.5 - 0.05; // for debug
        u = Vector3f::dot(p - d * _normal, ua);
        u = u / 4.43 - 0.5; // for debug

        h.set(t, material, get_normal(u, v, cos), material->get_color(u, v), p);
        return true;
    }

    Vector3f get_normal(float u, float v, float dec)
    {
        if (ua.squaredLength() < 1e-6)
        {
            if (dec < 0)
                return _normal;
            else
                return -_normal;
        }
        float res, u_grad, v_grad;

        material->update_normal(u, v, u_grad, v_grad, res);
        if (fabs(res) < 1e-6) //不涉及凹凸纹理贴图
        {
            if (dec < 0)
                return _normal;
            else
                return -_normal;
        }
        Vector3f new_normal = Vector3f::cross(ua + _normal * u_grad, Vector3f::UP + _normal * v_grad);
        new_normal.normalize();
        return new_normal;
    }

    // following code essential
    Vector3f min() const override
    {
        return -Vector3f(fabs(_normal.x()) < 1 - 1e-6,
                         fabs(_normal.y()) < 1 - 1e-6,
                         fabs(_normal.z()) < 1 - 1e-6) *
                   INF +
               _normal * d;
    }
    Vector3f max() const override
    {
        return Vector3f(fabs(_normal.x()) < 1 - 1e-6,
                        fabs(_normal.y()) < 1 - 1e-6,
                        fabs(_normal.z()) < 1 - 1e-6) *
                   INF +
               _normal * d;
    }

protected:
    Vector3f _normal;
    float d;
    Vector3f ua;
};

#endif // PLANE_H
