#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

#include "scene_parser.hpp"
using namespace std;
int main(int argc, char *argv[])
{
    for (int argNum = 1; argNum < argc; ++argNum)
    {
        std::cout << "Argument " << argNum << " is: " << argv[argNum]
                  << std::endl;
    }
    if (argc != 6)
    {
        std::cout << "Usage: ./bin/PA1 <inputfile> <output_dir> <iterations> <photons> <ckpt_interval>"
                  << endl;
        exit(1);
    }
    int iterations = atoi(argv[3]), photons = atoi(argv[4]), ckpt_interval = atoi(argv[5]);

    Scene scene(argv[1], argv[2], iterations, photons, ckpt_interval);
    scene.render();
    return 0;
}
