#pragma once

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>

#include <algorithm>
#include <string>
#include <assert.h>
#include <iostream>
using namespace std;

// TODO (PA2): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint
{
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve : public Object3D
{
protected:
    std::vector<Vector3f> controls;

public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {}

    bool intersect(const Ray &r, Hit &h) override
    {
        return false;
    }

    std::vector<Vector3f> &getControls()
    {
        return controls;
    }

    virtual void discretize(int resolution, std::vector<CurvePoint> &data) = 0;
};

class BernstainParam
{
public:
    BernstainParam(int _n, int _k, string _mode)
    {
        // attention, here n denotes the same meaning with basis.py, with n control points
        n = _n;
        k = _k;
        mode = _mode;
        knot = new float[n + k + 1];
        if (mode == BSPLINE_MODE)
        {
            for (int i = 0; i < n + k; i++)
                knot[i] = float(i) / float(n + k);
        }
        else if (mode == BEZIER_MODE)
        {
            assert(n == (k + 1));
            for (int i = 0; i < n; i++)
                knot[i] = 0.0;
            for (int i = n; i < 2 * n; i++)
                knot[i] = 1.0;
        }
        else
        {
            std::cout << "invalid mode" << std::endl;
            exit(-1);
        }
    }
    ~BernstainParam()
    {
        delete[] knot;
    }
    int get_valid_range(float *t_range, int res)
    {
        // std::cout << "calling get valid range" << std::endl;
        if (mode == BSPLINE_MODE)
        {
            float step = (1.0 / float(n + k)) / float(res);
            int cnt = 0;
            for (int i = k; i < n; i++)
            {
                for (int j = 0; j < res; j++)
                    t_range[cnt++] = knot[i] + j * step;
            }
            t_range[cnt] = knot[n];
            assert(cnt == (n - k) * res);
            return (n - k) * res + 1;
        }
        else
        {
            for (int i = 0; i < res; i++)
                t_range[i] = float(i) / float(res);
            t_range[res] = 1.0;
            return res + 1;
        }
    }
    int get_coeff(float *pt, float *dpt, float mu)
    {
        // std::cout << "calling get coeff" << std::endl;
        int bpos = get_bpos(mu);
        for (int i = 0; i < k + 2; i++)
            pt[i] = 0.0;
        for (int i = 0; i < k + 1; i++)
            dpt[i] = 1.0;
        pt[k] = 1.0;
        for (int p = 1; p < k + 1; p++)
        {
            for (int ii = k - p; ii < k + 1; ii++)
            {
                float w1, dw1;
                float w2, dw2;
                int i = ii + bpos - k;
                if (get_knot_pad(i + p) == get_knot_pad(i))
                {
                    w1 = mu;
                    dw1 = 1.0;
                }
                else
                {
                    w1 = (mu - get_knot_pad(i)) / (get_knot_pad(i + p) - get_knot_pad(i));
                    dw1 = 1.0 / (get_knot_pad(i + p) - get_knot_pad(i));
                }
                if (get_knot_pad(i + p + 1) == get_knot_pad(i + 1))
                {
                    w2 = 1 - mu;
                    dw2 = -1.0;
                }
                else
                {
                    w2 = (get_knot_pad(i + p + 1) - mu) / (get_knot_pad(i + p + 1) - get_knot_pad(i + 1));
                    dw2 = -1.0 / (get_knot_pad(i + p + 1) - get_knot_pad(i + 1));
                }
                if (p == k)
                    dpt[ii] = (dw1 * pt[ii] + dw2 * pt[ii + 1]) * p;
                pt[ii] = w1 * pt[ii] + w2 * pt[ii + 1];
            }
        }

        int lsk = bpos - k;
        int rsk = n - bpos - 1;
        if (lsk < 0 || rsk < 0)
        {
            std::cout << "lsk or rsk error" << std::endl;
            exit(-1);
        }
        return lsk;
    }

private:
    int n;
    int k;
    float *knot;
    string mode;
    static string BSPLINE_MODE;
    static string BEZIER_MODE;
    float get_knot_pad(int rank);
    int get_bpos(float mu)
    {
        // std::cout << "calling get bpos" << std::endl;
        assert((mu >= knot[0]) && (mu <= knot[n + k]));
        int curr_max = -1;
        for (int i = 0; i < n; i++)
        {
            if (mu >= knot[i])
            {
                curr_max = i;
            }
            else
            {
                break;
            }
        }
        return curr_max;
    }
};

string BernstainParam::BSPLINE_MODE = string("BSPLINE_MODE");
string BernstainParam::BEZIER_MODE = string("BEZIER_MODE");
float BernstainParam::get_knot_pad(int rank)
{
    // std::cout << "calling get knot pad" << std::endl;
    if ((rank < (n + k + 1)) && (rank >= 0))
        return knot[rank];
    else if ((rank >= (n + k + 1)) && (rank < n + 2 * k + 1))
        return knot[n + k];
    else
    {
        std::cout << "Invalid pad rank " << rank << std::endl;
        exit(-1);
    }
}

class BezierCurve : public Curve
{
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points)
    {
        if (points.size() < 4 || points.size() % 3 != 1)
        {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override
    {
        data.clear();
        // TODO (PA2): fill in data vector
        string mode = "BEZIER_MODE";
        // std::cout << "discretize BEZIER" << std::endl;
        int n = controls.size();
        int k = n - 1;
        BernstainParam bp(n, n - 1, mode);
        float *t_range;
        t_range = new float[(n - k) * resolution + 1];
        int t_range_len = bp.get_valid_range(t_range, resolution);
        for (int i = 0; i < t_range_len; i++)
        {
            // std::cout << "inside dis " << i << " and all lenth is " << t_range_len << std::endl;
            float *pt, *dpt;
            pt = new float[k + 2];
            dpt = new float[k + 1];
            int lsk = bp.get_coeff(pt, dpt, t_range[i]);
            CurvePoint curr;
            curr.T = Vector3f::ZERO;
            curr.V = Vector3f::ZERO;
            for (int j = lsk; j <= lsk + n - 1; j++)
            {
                curr.V = curr.V + pt[j - lsk] * controls[j];
                curr.T = curr.T + dpt[j - lsk] * controls[j];
            }
            curr.T.normalize();
            data.push_back(curr);
            delete[] pt;
            delete[] dpt;
        }
        delete[] t_range;
    }

protected:
};

class BsplineCurve : public Curve
{
public:
    BsplineCurve(const std::vector<Vector3f> &points) : Curve(points)
    {
        if (points.size() < 4)
        {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override
    {
        data.clear();
        // TODO (PA2): fill in data vector
        string mode = "BSPLINE_MODE";
        // std::cout << "discretize BSPLINE" << std::endl;
        int n = controls.size();
        int k = 3;
        BernstainParam bp(n, 3, mode);
        float *t_range;
        t_range = new float[(n - k) * resolution + 1];
        int t_range_len = bp.get_valid_range(t_range, resolution);
        for (int i = 0; i < t_range_len; i++)
        {
            // std::cout << "inside dis " << i << " and all lenth is " << t_range_len << std::endl;
            float *pt, *dpt;
            pt = new float[k + 2];
            dpt = new float[k + 1];
            int lsk = bp.get_coeff(pt, dpt, t_range[i]);
            CurvePoint curr;
            curr.T = Vector3f::ZERO;
            curr.V = Vector3f::ZERO;
            for (int j = lsk; j <= lsk + 3; j++)
            {
                curr.V = curr.V + pt[j - lsk] * controls[j];
                curr.T = curr.T + dpt[j - lsk] * controls[j];
            }
            curr.T.normalize();
            data.push_back(curr);
            delete[] pt;
            delete[] dpt;
        }
        delete[] t_range;
    }

protected:
};
