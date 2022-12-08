#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <float.h>

#include <cmath>
#include "utils.hpp"
#include <vecmath.h>
class Camera
{
public:
    Camera(const Vector3f &center, const Vector3f &direction,
           const Vector3f &up, int imgW, int imgH)
    {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

class PerspectiveCamera : public Camera
{
public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle, float _focal_lenth = 20.0f, float _aperture = 1.0f) : Camera(center, direction, up, imgW, imgH)
    {
        // angle is in radian.
        focal_len = _focal_lenth;
        aperture = _aperture;
        radian_angle = angle;
        cx = imgW / 2.0;
        cy = imgH / 2.0;
        fx = float(imgH) / (2.0 * tan(radian_angle / 2.0));
        fy = fx;
    }

    Ray generateRay(const Vector2f &point) override //加上了光圈，实现景深效果
    {
        //先转化为实际大小（focal length真的对应于focal length）
        Vector3f dr_pre((point.x() - cx) * focal_len / fx + RND * aperture, (cy - point.y()) * focal_len / fy + RND * aperture, focal_len);
        dr_pre.normalize();
        Matrix3f rotation(horizontal, -up, direction);

        Ray r(center, rotation * dr_pre);
        return r;
    }

protected:
    // Perspective intrinsics
    float cx, cy;
    float fx, fy;
    float radian_angle, aperture, focal_len;
};
#endif // CAMERA_H
