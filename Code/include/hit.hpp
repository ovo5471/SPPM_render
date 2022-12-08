#ifndef HIT_H
#define HIT_H

#include <vecmath.h>
#include "ray.hpp"
class Material;

#ifndef INF
#define INF 0x7fffffff
#endif

#ifndef EPS_R2
#define EPS_R2 1E-4
#endif

class Hit
{
public:
    Hit()
    {
        material = nullptr;
        t = INF;
        r2 = EPS_R2;
        n = 0;
        normal = fluxLight = flux = color = Vector3f::ZERO;
        attenuation = Vector3f(1, 1, 1);
    }

    Hit(float _t, Material *m, const Vector3f &_n)
    {
        t = _t;
        material = m;
        normal = _n;
        r2 = 1E-4;
        attenuation = Vector3f(1);
        fluxLight = flux = color = Vector3f::ZERO;
        n = 0;
    }

    float getT() const { return t; }

    Material *getMaterial() const { return material; }

    const Vector3f &get_normalal() const { return normal; }

    void reset_t()
    {
        t = INF;
    }
    void set(float _t, Material *m, const Vector3f &n, const Vector3f &c, const Vector3f &_p)
    {
        t = _t;
        material = m;
        normal = n;
        color = c;
        p = _p;
    }

    float t, r2;
    Material *material;
    Vector3f normal, color, flux, fluxLight, attenuation;
    Vector3f p;
    int n;
};

#endif // HIT_H
