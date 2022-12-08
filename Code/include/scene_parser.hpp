#ifndef SCENE_H
#define SCENE_H

#include <vecmath.h>

#include <cassert>
#include <vector>

class Camera;
class Material;
class Object3D;
class Group;
class Sphere;
class Plane;
class Triangle;
class Transform;
class Mesh;
class Curve;
class RevSurface;
class Image;
class HitKDTree;
class Hit;
#define MAX_PARSER_TOKEN_LENGTH 1024
#define DIFUSSION_TYPE 1
#define SPECULATION_TYPE 2
#define REFRACTION_TYPE 3

using namespace std;

class Scene
{
public:
    Scene() = delete;
    Scene(const char *filename, const char *dir, int _iters, int _photons, int _ckpt_interval, int _max_depth = 20);

    ~Scene();

    Camera *getCamera() const { return camera; }

    Vector3f getBackgroundColor() const { return background_color; }

    int getNumMaterials() const { return num_materials; }

    Material *getMaterial(int i) const
    {
        assert(i >= 0 && i < num_materials);
        return materials[i];
    }

    Group *getGroup() const { return group; }

    void render();

private:
    void parseFile();
    void parsePerspectiveCamera();
    void parseBackground();
    void parseLights();
    void parseDirectionalLight();
    void parsePointLight();

    int MT_sample(float x, float y, float z);

    void get_hit_points(int width, int height, vector<Hit *> &hitPoints); //获取visible points
    void allocate_photons(int curr_iter, int photon_allocated, HitKDTree *hitKDTree);
    void save_img(Image *img, const char *save_pth, int ckpt);
    void parseMaterials();
    Material *parseMaterial();
    Object3D *parseObject(char token[MAX_PARSER_TOKEN_LENGTH]);
    Group *parseGroup();
    Sphere *parseSphere();
    Plane *parsePlane();
    Triangle *parseTriangle();
    Mesh *parseTriangleMesh();
    Transform *parseTransform();
    Curve *parseBezierCurve();
    Curve *parseBsplineCurve();
    RevSurface *parseRevSurface();

    int getToken(char token[MAX_PARSER_TOKEN_LENGTH]);

    Vector3f readVector3f();

    float readFloat();
    int readInt();

    FILE *file;
    Camera *camera;
    Vector3f background_color;
    int num_lights;
    vector<Object3D *> illuminants; //可以发光的物体
    int num_materials;
    Material **materials;
    Material *current_material;
    Group *group;
    std::vector<Object3D *> objs;

    const char *output_dir;
    int iters, photons, ckpt_interval;
    int max_depth;
};

#endif // SCENE_H
