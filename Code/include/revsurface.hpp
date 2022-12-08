#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP
// 参数曲面
#include <tuple>

#include "bound.hpp"
#include "curve.hpp"
#include "object3d.hpp"
#include "triangle.hpp"
#define RESOLUTION_CURVE 100
#define ROTATION_STEP 40
class RevSurface : public Object3D
{
    Curve *pCurve;
    AABB aabb;
    // Definition for drawable surface.
    typedef std::tuple<unsigned, unsigned, unsigned> Tup3u;
    std::vector<Triangle> triangles;
    // Surface is just a struct that contains vertices, normals, and
    // faces.  VV[i] is the position of vertex i, and VN[i] is the normal
    // of vertex i.  A face is a triple i,j,k corresponding to a triangle
    // with (vertex i, normal i), (vertex j, normal j), ...
    // Currently this struct is computed every time when canvas refreshes.
    // You can store this as member function to accelerate rendering.
public:
    RevSurface(Curve *pCurve, Material *material)
        : pCurve(pCurve), Object3D(material)
    {
        // Check flat.
        for (const auto &cp : pCurve->getControls())
        {
            if (cp.z() != 0.0)
            {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        get_mesh();
    }

    ~RevSurface() override { delete pCurve; }

    inline bool intersect(const Ray &r, Hit &h) override
    {
        // return meshIntersect(r, h, tmin);
        float t;
        if (!aabb.intersect(r, t, 0) || t > h.getT())
            return false;

        //对所有的triangle求一次交
        bool itsc = false;
        for (int i = 0; i < triangles.size(); i++)
            if (triangles[i].intersect(r, h))
                itsc = true;
        return itsc;
    }

    Vector3f min() const override { return aabb.get_min(); }
    Vector3f max() const override { return aabb.get_max(); }

private:
    void get_mesh() //在产生三角面片的同时还要产生aabb的包围盒，方便求交快速判定
    {
        // Definition for drawable surface.
        // Surface is just a struct that contains vertices, normals, and
        // faces.  VV[i] is the position of vertex i, and VN[i] is the normal
        // of vertex i.  A face is a triple i,j,k corresponding to a triangle
        // with (vertex i, normal i), (vertex j, normal j), ...
        // Currently this struct is computed every time when canvas refreshes.
        // You can store this as member function to accelerate rendering.

        struct Surface
        {
            std::vector<Vector3f> VV;
            std::vector<Vector3f> VN;
            std::vector<Tup3u> VF;
            std::vector<float> u_coord;
            std::vector<float> v_coord;
        } surface;

        std::vector<CurvePoint> curvePoints;
        pCurve->discretize(RESOLUTION_CURVE, curvePoints);
        const int steps = ROTATION_STEP;
        float curr_u, curr_v;
        for (unsigned int ci = 0; ci < curvePoints.size(); ++ci)
        {
            curr_v = float(curvePoints.size() - 1 - ci) / float(curvePoints.size() - 1);
            const CurvePoint &cp = curvePoints[ci];
            for (unsigned int i = 0; i < steps; ++i)
            {
                curr_u = float(i) / float(steps);
                surface.u_coord.push_back(curr_u);
                surface.v_coord.push_back(curr_v);
                float t = (float)i / steps;
                Quat4f rot;
                rot.setAxisAngle(t * 2 * 3.14159, Vector3f::UP);
                Vector3f pnew = Matrix3f::rotation(rot) * cp.V;
                Vector3f pNormal = Vector3f::cross(cp.T, -Vector3f::FORWARD);
                Vector3f nnew = Matrix3f::rotation(rot) * pNormal;
                surface.VV.push_back(pnew);
                surface.VN.push_back(nnew);
                aabb.update(pnew);
                int i1 = (i + 1 == steps) ? 0 : i + 1;
                if (ci != curvePoints.size() - 1)
                {
                    surface.VF.emplace_back((ci + 1) * steps + i, ci * steps + i1, ci * steps + i);
                    surface.VF.emplace_back((ci + 1) * steps + i, (ci + 1) * steps + i1, ci * steps + i1);
                }
            }
        }

        for (unsigned i = 0; i < surface.VF.size(); i++) //创建三角面片，注意一下V和N的对应关系
        {
            Triangle t(surface.VV[std::get<0>(surface.VF[i])], surface.VV[std::get<1>(surface.VF[i])], surface.VV[std::get<2>(surface.VF[i])], material);
            t.set_normals(surface.VN[std::get<0>(surface.VF[i])], surface.VN[std::get<1>(surface.VF[i])], surface.VN[std::get<2>(surface.VF[i])]);
            t.set_uvs(Vector2f(surface.u_coord[std::get<0>(surface.VF[i])], surface.v_coord[std::get<0>(surface.VF[i])]),
                      Vector2f(surface.u_coord[std::get<1>(surface.VF[i])], surface.v_coord[std::get<1>(surface.VF[i])]),
                      Vector2f(surface.u_coord[std::get<2>(surface.VF[i])], surface.v_coord[std::get<2>(surface.VF[i])])); //设置UV坐标
            triangles.push_back(t);
        }
    }
};

#endif // REVSURFACE_HPP
