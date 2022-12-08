#include "scene_parser.hpp"

#include <vecmath.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>

#include "camera.hpp"
#include "curve.hpp"
#include "group.hpp"
#include "material.hpp"
#include "mesh.hpp"
#include "object3d.hpp"
#include "plane.hpp"
#include "revsurface.hpp"
#include "sphere.hpp"
#include "transform.hpp"
#include "triangle.hpp"
#include "image.hpp"
#include "hit.hpp"
#include "hit_kdtree.hpp"

#ifndef INF
#define INF 0x7fffffff;
#endif
using namespace std;
#define DegreesToRadians(x) ((M_PI * x) / 180.0f)

Scene::Scene(const char *filename, const char *dir, int _iters, int _photons, int _ckpt_interval, int _max_depth)
{
    // init parameters for rendering
    output_dir = dir;
    iters = _iters;
    photons = _photons;
    ckpt_interval = _ckpt_interval;
    max_depth = _max_depth;

    // initialize some reasonable default values
    group = nullptr;
    camera = nullptr;
    background_color = Vector3f(0.5, 0.5, 0.5);
    num_lights = 0;
    num_materials = 0;
    materials = nullptr;
    current_material = nullptr;

    // parse the file
    assert(filename != nullptr);
    const char *ext = &filename[strlen(filename) - 4];

    if (strcmp(ext, ".txt") != 0)
    {
        printf("wrong file name extension\n");
        exit(0);
    }
    file = fopen(filename, "r");

    if (file == nullptr)
    {
        printf("cannot open scene file\n");
        exit(0);
    }
    parseFile();
    fclose(file);
    file = nullptr;

    illuminants = group->illuminants;
    assert(illuminants.size() != 0);

    if (num_lights == 0)
    {
        printf("WARNING:    No lights specified\n");
    }
}

Scene::~Scene()
{
    delete group;
    delete camera;

    int i;
    for (i = 0; i < num_materials; i++)
    {
        delete materials[i];
    }
    delete[] materials;
}

// ====================================================================
// ====================================================================

void Scene::parseFile()
{
    //
    // at the top level, the scene can have a camera,
    // background color and a group of objects
    // (we add lights and other things in future assignments)
    //
    char token[MAX_PARSER_TOKEN_LENGTH];
    while (getToken(token))
    {
        if (!strcmp(token, "PerspectiveCamera"))
        {
            parsePerspectiveCamera();
        }
        else if (!strcmp(token, "Background"))
        {
            parseBackground();
        }
        else if (!strcmp(token, "Lights"))
        {
            parseLights();
        }
        else if (!strcmp(token, "Materials"))
        {
            parseMaterials();
        }
        else if (!strcmp(token, "Group"))
        {
            group = parseGroup();
        }
        else
        {
            printf("Unknown token in parseFile: '%s'\n", token);
            exit(0);
        }
    }
}

// ====================================================================
// ====================================================================

void Scene::parsePerspectiveCamera()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    Vector3f center, direction, up;
    float angle, focalLength = 1, aperture = 0;
    int width, height;
    // read in the camera parameters
    getToken(token);
    assert(!strcmp(token, "{"));
    while (true)
    {
        getToken(token);
        if (!strcmp(token, "center"))
        {
            center = readVector3f();
        }
        else if (!strcmp(token, "direction"))
        {
            direction = readVector3f();
        }
        else if (!strcmp(token, "up"))
        {
            up = readVector3f();
        }
        else if (!strcmp(token, "angle"))
        {
            float angle_degrees = readFloat();
            angle = DegreesToRadians(angle_degrees);
        }
        else if (!strcmp(token, "width"))
        {
            width = readInt();
        }
        else if (!strcmp(token, "height"))
        {
            height = readInt();
        }
        else if (strcmp(token, "focalLength") == 0)
        {
            focalLength = readFloat();
        }
        else if (strcmp(token, "aperture") == 0)
        {
            aperture = readFloat();
        }
        else
        {
            assert(!strcmp(token, "}"));
            break;
        }
    }
    camera = new PerspectiveCamera(center, direction, up, width, height, angle,
                                   focalLength, aperture);
}

void Scene::parseBackground()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    // read in the background color
    getToken(token);
    assert(!strcmp(token, "{"));
    while (true)
    {
        getToken(token);
        if (!strcmp(token, "}"))
        {
            break;
        }
        else if (!strcmp(token, "color"))
        {
            background_color = readVector3f();
        }
        else
        {
            printf("Unknown token in parseBackground: '%s'\n", token);
            assert(0);
        }
    }
}

// ====================================================================
// ====================================================================

void Scene::parseLights()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    // read in the number of objects
    getToken(token);
    assert(!strcmp(token, "numLights"));
    num_lights = readInt();
    // read in the objects
    int count = 0;
    while (num_lights > count)
    {
        getToken(token);
        if (strcmp(token, "DirectionalLight") == 0)
        {
            parseDirectionalLight();
        }
        else if (strcmp(token, "PointLight") == 0)
        {
            parsePointLight();
        }
        else
        {
            printf("Unknown token in parseLight: '%s'\n", token);
            exit(0);
        }
        count++;
    }
    getToken(token);
    assert(!strcmp(token, "}"));
}

void Scene::parseDirectionalLight()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "direction"));
    Vector3f direction = readVector3f();
    getToken(token);
    assert(!strcmp(token, "color"));
    Vector3f color = readVector3f();
    getToken(token);
    assert(!strcmp(token, "}"));
}

void Scene::parsePointLight()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "position"));
    Vector3f position = readVector3f();
    getToken(token);
    assert(!strcmp(token, "color"));
    Vector3f color = readVector3f();
    getToken(token);
    assert(!strcmp(token, "}"));
}
// ====================================================================
// ====================================================================

void Scene::parseMaterials()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    // read in the number of objects
    getToken(token);
    assert(!strcmp(token, "numMaterials"));
    num_materials = readInt();
    materials = new Material *[num_materials];
    // read in the objects
    int count = 0;
    while (num_materials > count)
    {
        getToken(token);
        if (!strcmp(token, "Material"))
        {
            materials[count] = parseMaterial();
        }
        else
        {
            printf("Unknown token in parseMaterial: '%s'\n", token);
            exit(0);
        }
        count++;
    }
    getToken(token);
    assert(!strcmp(token, "}"));
}

Material *Scene::parseMaterial()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    char textureFile[MAX_PARSER_TOKEN_LENGTH] = {0},
         bumpFile[MAX_PARSER_TOKEN_LENGTH] = {0};
    Vector3f color(1, 1, 1), specularColor(0, 0, 0), emission(0, 0, 0),
        type(1, 0, 0);
    float shininess = 0, refr = 1;
    getToken(token);
    assert(!strcmp(token, "{"));
    while (true)
    {
        getToken(token);
        if (strcmp(token, "color") == 0 || strcmp(token, "diffuseColor") == 0)
        {
            color = readVector3f();
        }
        else if (strcmp(token, "specularColor") == 0)
        {
            specularColor = readVector3f();
        }
        else if (strcmp(token, "shininess") == 0)
        {
            shininess = readFloat();
        }
        else if (strcmp(token, "refr") == 0)
        {
            refr = readFloat();
        }
        else if (strcmp(token, "texture") == 0)
        {
            // Optional: read in texture and draw it.
            getToken(textureFile);
        }
        else if (strcmp(token, "bump") == 0)
        {
            // Optional: read in texture and draw it.
            getToken(bumpFile);
        }
        else if (strcmp(token, "type") == 0)
        {
            type = readVector3f();
            type = type / (type.x() + type.y() + type.z());
        }
        else if (strcmp(token, "emission") == 0)
        {
            emission = readVector3f();
        }
        else
        {
            assert(!strcmp(token, "}"));
            break;
        }
    }
    auto *answer = new Material(color, specularColor, shininess, emission, refr,
                                type, textureFile, bumpFile);
    return answer;
}

// ====================================================================
// ====================================================================

Object3D *Scene::parseObject(char token[MAX_PARSER_TOKEN_LENGTH])
{
    Object3D *answer = nullptr;
    if (!strcmp(token, "Group"))
    {
        answer = (Object3D *)parseGroup();
    }
    else if (!strcmp(token, "Sphere"))
    {
        answer = (Object3D *)parseSphere();
    }
    else if (!strcmp(token, "Plane"))
    {
        answer = (Object3D *)parsePlane();
    }
    else if (!strcmp(token, "Triangle"))
    {
        answer = (Object3D *)parseTriangle();
    }
    else if (!strcmp(token, "TriangleMesh"))
    {
        answer = (Object3D *)parseTriangleMesh();
    }
    else if (!strcmp(token, "Transform"))
    {
        answer = (Object3D *)parseTransform();
    }
    else if (!strcmp(token, "BezierCurve"))
    {
        answer = (Object3D *)parseBezierCurve();
    }
    else if (!strcmp(token, "BsplineCurve"))
    {
        answer = (Object3D *)parseBsplineCurve();
    }
    else if (!strcmp(token, "RevSurface"))
    {
        answer = (Object3D *)parseRevSurface();
    }
    else
    {
        printf("Unknown token in parseObject: '%s'\n", token);
        exit(0);
    }
    return answer;
}

// ====================================================================
// ====================================================================

Group *Scene::parseGroup()
{
    //
    // each group starts with an integer that specifies
    // the number of objects in the group
    //
    // the material index sets the material of all objects which follow,
    // until the next material index (scoping for the materials is very
    // simple, and essentially ignores any tree hierarchy)
    //
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));

    // read in the number of objects
    getToken(token);
    assert(!strcmp(token, "numObjects"));
    int num_objects = readInt();
    // read in the objects
    int count = 0;
    while (num_objects > count)
    {
        getToken(token);
        if (!strcmp(token, "MaterialIndex"))
        {
            // change the current material
            int index = readInt();
            assert(index >= 0 && index <= getNumMaterials());
            current_material = getMaterial(index);
        }
        else
        {
            Object3D *object = parseObject(token);
            assert(object != nullptr);
            objs.push_back(object);
            count++;
        }
    }
    getToken(token);
    assert(!strcmp(token, "}"));

    // return the group
    return new Group(objs);
}

// ====================================================================
// ====================================================================

Sphere *Scene::parseSphere()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "center"));
    Vector3f center = readVector3f();
    getToken(token);
    assert(!strcmp(token, "radius"));
    float radius = readFloat();
    getToken(token);
    assert(!strcmp(token, "}"));
    assert(current_material != nullptr);
    return new Sphere(center, radius, current_material);
}

Plane *Scene::parsePlane()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "normal"));
    Vector3f normal = readVector3f();
    getToken(token);
    assert(!strcmp(token, "offset"));
    float offset = readFloat();
    getToken(token);
    assert(!strcmp(token, "}"));
    assert(current_material != nullptr);
    return new Plane(normal, offset, current_material);
}

Triangle *Scene::parseTriangle()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "vertex0"));
    Vector3f v0 = readVector3f();
    getToken(token);
    assert(!strcmp(token, "vertex1"));
    Vector3f v1 = readVector3f();
    getToken(token);
    assert(!strcmp(token, "vertex2"));
    Vector3f v2 = readVector3f();
    getToken(token);
    assert(!strcmp(token, "}"));
    assert(current_material != nullptr);
    return new Triangle(v0, v1, v2, current_material);
}

Mesh *Scene::parseTriangleMesh()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    char filename[MAX_PARSER_TOKEN_LENGTH];
    // get the filename
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "obj_file"));
    getToken(filename);
    getToken(token);
    assert(!strcmp(token, "}"));
    const char *ext = &filename[strlen(filename) - 4];
    assert(!strcmp(ext, ".obj"));
    Mesh *answer = new Mesh(filename, current_material);

    return answer;
}

Curve *Scene::parseBezierCurve()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "controls"));
    vector<Vector3f> controls;
    while (true)
    {
        getToken(token);
        if (!strcmp(token, "["))
        {
            controls.push_back(readVector3f());
            getToken(token);
            assert(!strcmp(token, "]"));
        }
        else if (!strcmp(token, "}"))
        {
            break;
        }
        else
        {
            printf("Incorrect format for BezierCurve!\n");
            exit(0);
        }
    }
    Curve *answer = new BezierCurve(controls);
    return answer;
}

Curve *Scene::parseBsplineCurve()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "controls"));
    vector<Vector3f> controls;
    while (true)
    {
        getToken(token);
        if (!strcmp(token, "["))
        {
            controls.push_back(readVector3f());
            getToken(token);
            assert(!strcmp(token, "]"));
        }
        else if (!strcmp(token, "}"))
        {
            break;
        }
        else
        {
            printf("Incorrect format for BsplineCurve!\n");
            exit(0);
        }
    }
    Curve *answer = new BsplineCurve(controls);
    return answer;
}

RevSurface *Scene::parseRevSurface()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    getToken(token);
    assert(!strcmp(token, "{"));
    getToken(token);
    assert(!strcmp(token, "profile"));
    Curve *profile;
    getToken(token);
    if (!strcmp(token, "BezierCurve"))
    {
        profile = parseBezierCurve();
    }
    else if (!strcmp(token, "BsplineCurve"))
    {
        profile = parseBsplineCurve();
    }
    else
    {
        printf("Unknown profile type in parseRevSurface: '%s'\n", token);
        exit(0);
    }
    getToken(token);
    assert(!strcmp(token, "}"));
    auto *answer = new RevSurface(profile, current_material);
    return answer;
}

Transform *Scene::parseTransform()
{
    char token[MAX_PARSER_TOKEN_LENGTH];
    Matrix4f matrix = Matrix4f::identity();
    Object3D *object = nullptr;
    getToken(token);
    assert(!strcmp(token, "{"));
    // read in transformations:
    // apply to the LEFT side of the current matrix (so the first
    // transform in the list is the last applied to the object)
    getToken(token);

    while (true)
    {
        if (!strcmp(token, "Scale"))
        {
            Vector3f s = readVector3f();
            matrix = matrix * Matrix4f::scaling(s[0], s[1], s[2]);
        }
        else if (!strcmp(token, "UniformScale"))
        {
            float s = readFloat();
            matrix = matrix * Matrix4f::uniformScaling(s);
        }
        else if (!strcmp(token, "Translate"))
        {
            matrix = matrix * Matrix4f::translation(readVector3f());
        }
        else if (!strcmp(token, "XRotate"))
        {
            matrix = matrix * Matrix4f::rotateX(DegreesToRadians(readFloat()));
        }
        else if (!strcmp(token, "YRotate"))
        {
            matrix = matrix * Matrix4f::rotateY(DegreesToRadians(readFloat()));
        }
        else if (!strcmp(token, "ZRotate"))
        {
            matrix = matrix * Matrix4f::rotateZ(DegreesToRadians(readFloat()));
        }
        else if (!strcmp(token, "Rotate"))
        {
            getToken(token);
            assert(!strcmp(token, "{"));
            Vector3f axis = readVector3f();
            float degrees = readFloat();
            float radians = DegreesToRadians(degrees);
            matrix = matrix * Matrix4f::rotation(axis, radians);
            getToken(token);
            assert(!strcmp(token, "}"));
        }
        else if (!strcmp(token, "Matrix4f"))
        {
            Matrix4f matrix2 = Matrix4f::identity();
            getToken(token);
            assert(!strcmp(token, "{"));
            for (int j = 0; j < 4; j++)
            {
                for (int i = 0; i < 4; i++)
                {
                    float v = readFloat();
                    matrix2(i, j) = v;
                }
            }
            getToken(token);
            assert(!strcmp(token, "}"));
            matrix = matrix2 * matrix;
        }
        else
        {
            // otherwise this must be an object,
            // and there are no more transformations
            object = parseObject(token);
            break;
        }
        getToken(token);
    }

    assert(object != nullptr);
    getToken(token);
    assert(!strcmp(token, "}"));
    return new Transform(matrix, object);
}

// ====================================================================
// ====================================================================

int Scene::getToken(char token[MAX_PARSER_TOKEN_LENGTH])
{
    // for simplicity, tokens must be separated by whitespace
    assert(file != nullptr);
    int success = fscanf(file, "%s ", token);
    if (success == EOF)
    {
        token[0] = '\0';
        return 0;
    }
    return 1;
}

Vector3f Scene::readVector3f()
{
    float x, y, z;
    int count = fscanf(file, "%f %f %f", &x, &y, &z);
    if (count != 3)
    {
        printf("Error trying to read 3 floats to make a Vector3f\n");
        assert(0);
    }
    return Vector3f(x, y, z);
}

float Scene::readFloat()
{
    float answer;
    int count = fscanf(file, "%f", &answer);
    if (count != 1)
    {
        printf("Error trying to read 1 float\n");
        assert(0);
    }
    return answer;
}

int Scene::readInt()
{
    int answer;
    int count = fscanf(file, "%d", &answer);
    if (count != 1)
    {
        printf("Error trying to read 1 int\n");
        assert(0);
    }
    return answer;
}

// ====================================================================
// ====================================================================

int Scene::MT_sample(float x, float y, float z)
{
    float s = RND2;
    if (s < x) // diffusion
        return DIFUSSION_TYPE;
    else if (s < x + y) // speculation
        return SPECULATION_TYPE;
    else
        return REFRACTION_TYPE; // refraction
}

// ====================================================================
// ====================================================================
void Scene::get_hit_points(int width, int height, vector<Hit *> &hitPoints)
{
#pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            int depth = 0;
            Ray camRay = camera->generateRay(Vector2f(x + RND, y + RND)); //获取相机视线，抗锯齿
            Hit *hitcoord = hitPoints[x * height + y];
            hitcoord->reset_t();
            Vector3f attenuation(1, 1, 1);
            while ((depth < max_depth) && (attenuation.max() > 1e-3)) //退出循环条件：深度到达最大次数，或者衰减达到预设成都
            {
                depth++;
                hitcoord->t = INF;                        //这里必须设置为INF,因为本质上反射和折射相当于移动相机位置，需要重新计算相交
                if (!group->intersect(camRay, *hitcoord)) //如果没有相交，设置为背景颜色乘以衰减
                {
                    hitcoord->fluxLight += hitcoord->attenuation * background_color;
                    break;
                }
                Material *material = (*hitcoord).material;
                int random_ray_type = MT_sample(material->type.x(), material->type.y(), material->type.z());
                if (random_ray_type == DIFUSSION_TYPE)
                { // Diffuse
                    hitcoord->attenuation = attenuation * hitcoord->color;
                    hitcoord->fluxLight += hitcoord->attenuation * material->emission;
                    break;
                }
                else if (random_ray_type == SPECULATION_TYPE) // speculation
                {
                    camRay.origin += camRay.direction * hitcoord->getT();
                    camRay.direction = (camRay.direction + (-2) * Vector3f::dot(camRay.direction, (hitcoord->normal)) * (hitcoord->normal));
                    camRay.direction.normalize();
                }
                else
                {
                    //折射光线近似求解借鉴https://github.com/Guangxuan-Xiao/THU-Computer-Graphics-2020/tree/master/project
                    camRay.origin += camRay.direction * hitcoord->getT();
                    float n = material->refr;
                    float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                    if (Vector3f::dot((hitcoord->normal), camRay.direction) > 0)
                    { // inside the medium
                        (hitcoord->normal).negate();
                        n = 1 / n;
                    }
                    n = 1 / n;
                    float cost1 =
                        -Vector3f::dot((hitcoord->normal), camRay.direction); // cosine theta_1
                    float cost2 =
                        1.0 - n * n * (1.0 - cost1 * cost1); // cosine theta_2
                    float Rprob =
                        R0 + (1.0 - R0) * pow(1.0 - cost1,
                                              5.0); // Schlick-approximation
                    if (cost2 > 0 && RND2 > Rprob)
                    { // refraction direction
                        camRay.direction = ((camRay.direction * n) + ((hitcoord->normal) * (n * cost1 - sqrt(cost2)))).normalized();
                    }
                    else
                    { // speculation direction
                        camRay.direction = (camRay.direction + (hitcoord->normal) * (cost1 * 2)).normalized();
                    }
                }
                attenuation = attenuation * hitcoord->color;
            }
        }
    }
}

// ====================================================================
// ====================================================================
void Scene::allocate_photons(int curr_iter, int photon_allocated, HitKDTree *hitKDTree)
{
#pragma omp parallel for schedule(dynamic, 1)
    for (int photon_rank = 0; photon_rank < photon_allocated; photon_rank++)
    {
        for (int illuminant_rank = 0; illuminant_rank < illuminants.size(); illuminant_rank++)
        {
            int depth = 0;
            Ray ray = illuminants[illuminant_rank]->randomRay(-1, (long long)curr_iter * photons + photon_allocated * illuminant_rank + photon_rank);
            Vector3f attenuation = illuminants[illuminant_rank]->material->emission * Vector3f(250, 250, 250);
            while ((++depth < max_depth) && (attenuation.max() > 1e-3))
            {
                Hit hit;
                if (!group->intersect(ray, hit))
                    break;
                Material *material = hit.material;
                int ray_pattern = MT_sample(material->type.x(), material->type.y(), material->type.z());
                ray.origin += ray.direction * hit.getT();
                if (ray_pattern == DIFUSSION_TYPE)
                {                                                                          // Diffuse
                    hitKDTree->update(hitKDTree->root, hit.p, attenuation, ray.direction); //完成一次送光
                    ray.direction = diffDir(hit.normal, -1, (long long)curr_iter * photons + photon_allocated * illuminant_rank + photon_rank + 13322);
                }
                else if (ray_pattern == SPECULATION_TYPE)
                {
                    ray.direction = (ray.direction + (-2) * Vector3f::dot(ray.direction, hit.normal) * hit.normal);
                    ray.direction.normalize();
                }
                else
                {
                    Vector3f N(hit.normal);
                    float n = material->refr;
                    float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
                    if (Vector3f::dot(N, ray.direction) > 0)
                    { // inside the medium
                        N.negate();
                        n = 1 / n;
                    }
                    n = 1 / n;
                    float cost1 =
                        -Vector3f::dot(N, ray.direction); // cosine theta_1
                    float cost2 =
                        1.0 - n * n * (1.0 - cost1 * cost1); // cosine theta_2
                    float Rprob =
                        R0 + (1.0 - R0) * pow(1.0 - cost1,
                                              5.0); // Schlick-approximation
                    if (cost2 > 0 && RND2 > Rprob)
                    { // refraction direction
                        ray.direction =
                            ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2))))
                                .normalized();
                    }
                    else
                    { // reflection direction
                        ray.direction =
                            (ray.direction + N * (cost1 * 2)).normalized();
                    }
                }
                attenuation = attenuation * hit.color;
            }
        }
    }
}

// ====================================================================
// ====================================================================
void Scene::render()
{
    //新建一张图片，先调用forward，然后调用backward函数赋予每一个光点颜色
    int width = camera->getWidth(), height = camera->getHeight();
    Image img(width, height);
    //确认创建画布，确认画布大小
    std::cout << "Canvas created.. Width=" << width << " Height=" << height << endl;

    //新建hit tree
    time_t start = time(NULL);
    vector<Hit *> hitpoints;
    for (int i = 0; i < width * height; i++)
        hitpoints.push_back(new Hit());
    HitKDTree *hitKDTree = nullptr;
    for (int curr_iter = 0; curr_iter < iters; curr_iter++)
    {
        //确认时间
        time_t curr_time = time(NULL);
        float interval = curr_time - start;
        float progress = float(curr_iter) / float(iters);
        std::cout << "Current iter=" << curr_iter + 1 << " Whole iter=" << iters << " time spent=" << interval << " progress=" << progress << endl;

        get_hit_points(camera->getWidth(), camera->getHeight(), hitpoints);

        if (hitKDTree != nullptr)
            delete hitKDTree;
        hitKDTree = new HitKDTree(&hitpoints);

        //每一个光源可以分配的光子数目
        int photon_allocated = photons / illuminants.size();
        allocate_photons(curr_iter, photon_allocated, hitKDTree);

        if ((curr_iter + 1) % ckpt_interval == 0) //到了checkpoint
        {
            for (int x = 0; x < width; ++x) //计算每一个像素点的光通量
                for (int y = 0; y < height; ++y)
                {
                    Hit *hit = hitpoints[x * height + y];
                    img.SetPixel(x, y, hit->flux / (M_PI * hit->r2 * photons * (curr_iter + 1)) + hit->fluxLight / (curr_iter + 1));
                }
            save_img(&img, output_dir, curr_iter);
        }
    }

    if (hitKDTree)
        delete hitKDTree;
    save_img(&img, output_dir, -1);
}

// ====================================================================
// ====================================================================

void Scene::save_img(Image *img, const char *save_pth, int ckpt)
{
    if (ckpt == -1)
        img->SaveBMP("result.bmp");
    else
    {
        char filename[50];
        sprintf(filename, "%s/ckpt-%d.bmp", save_pth, ckpt + 1);
        img->SaveBMP(filename);
    }
}