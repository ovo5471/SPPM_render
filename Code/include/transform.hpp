#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vecmath.h>
#include "object3d.hpp"
#include "bound.hpp"

// transforms a 3D point using a matrix, returning a 3D point
static Vector3f transformPoint(const Matrix4f &mat, const Vector3f &point)
{
    return (mat * Vector4f(point, 1)).xyz();
}

// transform a 3D directino using a matrix, returning a direction
static Vector3f transformDirection(const Matrix4f &mat, const Vector3f &dir)
{
    return (mat * Vector4f(dir, 0)).xyz();
}

class Transform : public Object3D
{
public:
    Transform() {}

    Transform(const Matrix4f &m, Object3D *obj) : o(obj), Object3D(obj->material)
    {
        transform = m.inverse();
        Vector3f prev_min(o->min());
        Vector3f prev_max(o->max());
        aabb.set(o->min(), o->max()); //不能直接对min和max进行旋转
    }

    void update_boundary(const Matrix4f &m)
    {
        //把aabb的八个顶点全部求出来，然后都对aabb进行update
        Vector3f *verts = aabb.get_8_vertices();
        aabb.init();
        for (int i = 0; i < 8; i++)
            aabb.update(transformPoint(m, verts[i]));
        delete[] verts;
    }

    ~Transform() {}

    bool intersect(const Ray &r, Hit &h) override
    {
        Vector3f trSource = transformPoint(transform, r.getOrigin());
        Vector3f trDirection = transformDirection(transform, r.getDirection());
        Ray tr(trSource, trDirection);
        if (!o->intersect(tr, h))
            return false;
        h.set(h.t, h.material,
              transformDirection(transform.transposed(), h.get_normalal())
                  .normalized(),
              h.color, h.t * r.direction + r.origin);
        return true;
    }

    Vector3f min() const override { return aabb.get_min(); }
    Vector3f max() const override { return aabb.get_max(); }

protected:
    Object3D *o; // un-transformed object
    Matrix4f transform;
    AABB aabb;
};

#endif // TRANSFORM_H
