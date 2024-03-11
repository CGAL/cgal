/**
 * @file   ui/vecmath.cpp
 * @author Gernot Walzl
 * @date   2013-06-15
 */

#include "ui/vecmath.h"

namespace ui {

void copy(const vec3f vec, vec3f& out) {
    for (unsigned int i = 0; i < 3; i++) {
        out[i] = vec[i];
    }
}

float length(const vec3f vec) {
    float result = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        result += vec[i] * vec[i];
    }
    result = sqrtf(result);
    return result;
}

void normalize(vec3f& vec) {
    float len = length(vec);
    for (unsigned int i = 0; i < 3; i++) {
        vec[i] /= len;
    }
}

float distance(const vec3f v1, const vec3f v2) {
    float result = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        float diff = v2[i] - v1[i];
        result += diff * diff;
    }
    result = sqrtf(result);
    return result;
}

void scale(vec3f& vec, float s) {
    for (unsigned int i = 0; i < 3; i++) {
        vec[i] *= s;
    }
}

float scalar(const vec3f v1, const vec3f v2) {
    float result = 0.0f;
    for (unsigned int i = 0; i < 3; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

void cross(const vec3f v1, const vec3f v2, vec3f& out) {
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

float angle(const vec3f v1, const vec3f v2) {
    float result = 0.0;
    vec3f n1;
    vec3f n2;
    copy(v1, n1);
    copy(v2, n2);
    normalize(n1);
    normalize(n2);
    float arg = 0.0f;
    for (unsigned int i = 0; i < 3; i++) {
        arg += n1[i] * n2[i];
    }
    if (arg <= -1.0f) {
        result = M_PI;
    } else if (arg >= 1.0f) {
        result = 0.0f;
    } else {
        result = acosf(arg);
    }
    return result;
}

}
