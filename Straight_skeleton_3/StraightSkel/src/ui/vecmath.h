/**
 * @file   ui/vecmath.h
 * @author Gernot Walzl
 * @date   2013-06-15
 */

#ifndef UI_VECMATH_H
#define UI_VECMATH_H

#include <cmath>
#include "ui/typedefs.h"

namespace ui {

void copy(const vec3f vec, vec3f& out);
float length(const vec3f vec);
void normalize(vec3f& vec);
float distance(const vec3f v1, const vec3f v2);
void scale(vec3f& vec, float s);
float scalar(const vec3f v1, const vec3f v2);
void cross(const vec3f v1, const vec3f v2, vec3f& out);
float angle(const vec3f v1, const vec3f v2);

}

#endif /* UI_VECMATH_H */
