#ifndef MODEL_H
#define MODEL_H

#include "hitable.h"
#include "vector.h"

#include <vector>

struct Triangle {
    Vec3f positions[3];
    Vec3f normals[3];
};

struct Model {
   std::vector<Triangle> triangles;
};

#endif // MODEL_H
