#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vecmath.h>

#include "object3d.hpp"
#include "utils.hpp"
#include "bound.hpp"

#include <cfloat>
#include <cmath>
#include <iostream>

using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle : public Object3D
{
public:
    Triangle() = delete;

    // a b c are three vertex positions of the triangle
    // triangle应该也可以通过aabb加速
    Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c,
             Material *m)
        : Object3D(m)
    {
        verts[0] = a;
        verts[1] = b;
        verts[2] = c;
        update_aabb();

        normal_set = false;
        uv_set = false;

        for (int i = 0; i < 3; i++)
            normals[i] = Vector3f::ZERO;

        normal = Vector3f::cross((verts[1] - verts[0]), (verts[2] - verts[0]));
        normal.normalize();
    }

    void update_aabb()
    {
        for (int i = 0; i < 3; i++)
            aabb.update(verts[i]);
    }

    bool intersect(const Ray &r, Hit &h) override
    {
        Vector3f r_origin = r.getOrigin();
        Vector3f r_dir = r.getDirection();

        Vector3f E1, E2, S;
        E1 = verts[0] - verts[1];
        E2 = verts[0] - verts[2];
        S = verts[0] - r_origin;
        Matrix3f M3(r_dir, E1, S);
        Matrix3f M4(r_dir, E1, E2);

        float m4_det = M4.determinant();
        float m3_det = M3.determinant();
        float gamma = m3_det / m4_det;
        if ((gamma < 0) || (gamma > 1))
            return false;

        Matrix3f M2(r_dir, S, E2);
        float m2_det = M2.determinant();
        float beta = m2_det / m4_det;
        if ((beta < 0) || (beta > 1) || (gamma + beta > 1))
            return false;
        float u = gamma, v = beta;

        Matrix3f M1(S, E1, E2, true);
        float m1_det = M1.determinant();
        float curr_t = m1_det / m4_det;

        if ((curr_t >= 0) && (curr_t <= h.getT()))
        {
            Vector3f p(r_origin + curr_t * r_dir);
            get_uv_coord(p, u, v);
            h.set(curr_t, material, get_normal(p, r_dir), material->get_color(u, v), p);
            return true;
        }
        return false;
    }

    void set_normals(const Vector3f &anorm, const Vector3f &bnorm,
                     const Vector3f &cnorm)
    {
        normal_set = true;
        normals[0] = anorm;
        normals[1] = bnorm;
        normals[2] = cnorm;
    }

    void set_uvs(const Vector2f &_at, const Vector2f &_bt, const Vector2f &_ct)
    {
        uv_set = true;
        uvs[0] = _at;
        uvs[1] = _bt;
        uvs[2] = _ct;
    }

    Vector3f min() const override { return aabb.get_min(); }
    Vector3f max() const override { return aabb.get_max(); }

private:
    Vector3f normal;
    Vector3f verts[3];
    Vector3f normals[3];
    Vector2f uvs[3];
    AABB aabb;
    bool normal_set = false;
    bool uv_set = false;
    void get_coe(float *weights, const Vector3f &p)
    {
        Vector3f edges[3];
        for (int i = 0; i < 3; i++)
            edges[i] = verts[i] - p;
        for (int i = 0; i < 3; i++)
            weights[i] = Vector3f::cross(edges[(i + 1) % 3], edges[(i + 2) % 3]).length();
    }

protected:
    Vector3f get_normal(const Vector3f &p, const Vector3f &r_dir)
    {
        if (!normal_set)
        {
            if (Vector3f::dot(r_dir, normal) > 0)
                return -normal;
            else
                return normal;
        }
        float weights[3];
        get_coe(weights, p);
        Vector3f res = Vector3f::ZERO;
        for (int i = 0; i < 3; i++)
            res = res + normals[i] * weights[i];
        res.normalize();

        return res;
    }

    void get_uv_coord(const Vector3f &p, float &u, float &v)
    {
        if (!uv_set)
            return;

        float weights[3];
        get_coe(weights, p);
        Vector2f res = Vector2f::ZERO;
        for (int i = 0; i < 3; i++)
            res = res + uvs[i] * weights[i] / (weights[0] + weights[1] + weights[2]);
        // Vector2f uv = (ra * uvs[0] + rb * uvs[1] + rc * uvs[2]) / (ra + rb + rc);
        u = res.x();
        v = res.y();
    }

    Ray randomRay(int axis = -1, long long int seed = 0) const override
    {
        return Ray(Vector3f::ZERO, Vector3f::ZERO);
    }
};

#endif // TRIANGLE_H
