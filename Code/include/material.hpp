#ifndef MATERIAL_H
#define MATERIAL_H
#include <cassert>
#include <assert.h>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <string>

#define BUMP_CHANNEL 1

using namespace std;

class Texture
{ // 纹理
public:
    unsigned char *pic;
    int w, h;
    Texture(const char *textureFile);
    Vector3f get_color(float u, float v) const;
    void get_corner(int u, int v, Vector3f &A, Vector3f &B, Vector3f &C, Vector3f &D) const;
    inline int clamp(int num, int lim)
    {
        return num % lim;
    }
};

Texture::Texture(const char *textureFile)
{
    pic = nullptr;
    if (strlen(textureFile) > 0)
    {
        int c;
        pic = stbi_load(textureFile, &w, &h, &c, 0);
        assert(c == 3);
    }
}

void Texture::get_corner(int u, int v, Vector3f &A, Vector3f &B, Vector3f &C, Vector3f &D) const
{
    int rank = u * 3 + v * w * 3;
    A = Vector3f(pic[rank + 0], pic[rank + 1], pic[rank + 2]) / 255.;
    if (u == w - 1)
        B = A;
    else
        B = Vector3f(pic[rank + 0 + 3], pic[rank + 1 + 3], pic[rank + 2 + 3]) / 255.;
    if (v == h - 1)
        C = A;
    else
        C = Vector3f(pic[rank + 0 + 3 + 3 * w], pic[rank + 1 + 3 + 3 * w], pic[rank + 2 + 3 + 3 * w]) / 255.;
    if (u == w - 1 || v == h - 1)
        D = A;
    else
        D = Vector3f(pic[rank + 0 + 3 * w], pic[rank + 1 + 3 * w], pic[rank + 2 + 3 * w]) / 255.;
}

Vector3f Texture::get_color(float u, float v) const
{
    if (!pic)
        return Vector3f::ZERO;
    u = float(w) * (u - floor(u));
    v = float(h) * (ceil(v) - v);
    int iu = (int)u, iv = (int)v;
    float alpha = u - iu, beta = v - iv;
    Vector3f ret, A, B, C, D;
    get_corner(iu, iv, A, B, C, D);
    ret = (1 - alpha) * (1 - beta) * A + alpha * (1 - beta) * B + (1 - alpha) * beta * D + alpha * beta * C; //双线性插值
    return ret;
}

class Bump
{ // 凹凸纹理
public:
    unsigned char *pic;
    int w, h;
    Bump(const char *textureFile);
    void get_grad(float u, float v, float &u_grad, float &v_grad, float res) const;
    inline int rank(float u, float v) const;
    inline float getGray(int idx) const;
};

inline int Bump::rank(float u, float v) const
{
    return (int(v * h) % h) * w * BUMP_CHANNEL + (int(u * w) % w) * BUMP_CHANNEL;
}

float Bump::getGray(int idx) const { return (pic[idx] / 255. - 0.5) * 2; }

Bump::Bump(const char *textureFile)
{
    pic = nullptr;
    if (strlen(textureFile) > 0)
    {
        int c;
        pic = stbi_load(textureFile, &w, &h, &c, 0);
        assert(c == BUMP_CHANNEL);
    }
}

void Bump::get_grad(float u, float v, float &u_grad, float &v_grad, float res) const
{
    if (!pic)
    {
        res = 0.0;
        return;
    }
    res = (pic[rank(u, v)] / 255. - 0.5) * 2;
    if (res < 1e-6)
    {
        res = 0.0;
        return;
    }
    u_grad = w * (getGray(rank(u + 1.0 / w, v)) - getGray(rank(u - 1.0 / w, v))) / 2.0;
    v_grad = h * (getGray(rank(u, v + 1.0 / h)) - getGray(rank(u, v - 1.0 / h))) / 2.0;
}

class Material
{
public:
    explicit Material(const Vector3f &_color,
                      const Vector3f &s_color = Vector3f::ZERO, float s = 0,
                      const Vector3f &e_color = Vector3f::ZERO, float refraction = 1,
                      Vector3f t = Vector3f(1, 0, 0),
                      const char *textureFile = "", const char *bumpFile = "")
        : texture(textureFile),
          bump(bumpFile)
    {

        color = _color;
        shininess = s;
        emission = e_color;
        refr = refraction;
        type = t;
        material_cnt++;
        if ((Vector3f::vec_sum(type) > 1.001) || (Vector3f::vec_sum(type < 0.999)))
        {
            std::cout << "Warning, type sum not 1. Check your material " << material_cnt << endl;
        }

        // assert((Vector3f::vec_sum(type) < 1.001) && (Vector3f::vec_sum(type > 0.999))); //确认type的和为1
    }

    virtual ~Material() = default;

    Vector3f get_color(float u, float v) const
    {
        if (!texture.pic)
            return color;
        else
            return texture.get_color(u, v);
    }

    float update_normal(float u, float v, float &u_grad, float &v_grad, float &res)
    {
        bump.get_grad(u, v, u_grad, v_grad, res);
    }

    Vector3f color;    // 颜色
    Vector3f emission; // 发光系数
    float shininess;   // 高光指数
    float refr;        // 折射率
    Vector3f type;     // 种类
    Texture texture;   // 纹理
    Bump bump;         //凹凸纹理
    static int material_cnt;
};

int Material::material_cnt = 0;

#endif // MATERIAL_H
