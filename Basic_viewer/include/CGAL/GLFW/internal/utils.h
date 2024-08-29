#ifndef CGAL_UTILS_EIGEN_H
#define CGAL_UTILS_EIGEN_H

#include <cmath>
#include <chrono>


#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

std::chrono::high_resolution_clock::time_point start;
std::chrono::high_resolution_clock::time_point stop;

std::chrono::microseconds duration;

#define MEASURE_TIME(func) \
        start = std::chrono::high_resolution_clock::now(); \
        func; \
        stop = std::chrono::high_resolution_clock::now(); \
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); \
        std::cout << "Time taken by " << #func << ": " << (duration.count() / 1000.0) << "ms" << std::endl; \

using vec2i = Eigen::Vector2i;
using vec2f = Eigen::Vector2f;
using vec3f = Eigen::Vector3f;
using vec4f = Eigen::Vector4f;

using mat3f = Eigen::Matrix3f;
using mat4f = Eigen::Matrix4f;

using quatf = Eigen::Quaternionf;

namespace utils
{
  bool inside_bounding_box_2d(const vec2f& point, const vec2f& pmin, const vec2f& pmax)
  {
    return point.x() >= pmin.x() && point.x() <= pmax.x() && point.y() >= pmin.y() && point.y() <= pmax.y(); 
  } 

  vec2f normalized_coordinates(const vec2f& coords, float width, float height)
  {
    vec2f result; 
    result.x() = coords.x() / width * 2 - 1;
    result.y() = -1 * (coords.y() / height * 2 - 1);

    return result; 
  }

  bool equal_float(float a, float b, float epsilon=0.001f)
  {
    return fabs(a-b) < epsilon;  
  }

  bool equal_vec3f(vec3f u, vec3f v, float epsilon=0.001f)
  {
    return equal_float(u.x(), v.x(), epsilon) && equal_float(u.y(), v.y(), epsilon) && equal_float(u.z(), v.z(), epsilon);  
  }

  float radians(float x)
  {
    return M_PI / 180 * x;
  }

  vec2f radians(const vec2f &v)
  {
    vec2f result;
    result << radians(v.x()),
        radians(v.y());
    return result;
  }

  float degrees(float x)
  {
    return 180 / M_PI * x;
  }

  mat4f euler_angle_xy(float const &angleX, float const &angleY)
  {
    float cosX = std::cos(angleX);
    float sinX = std::sin(angleX);
    float cosY = std::cos(angleY);
    float sinY = std::sin(angleY);

    mat4f result = mat4f::Identity();
    result(0, 0) = cosY;
    result(1, 0) = -sinX * -sinY;
    result(2, 0) = cosX * -sinY;
    result(1, 1) = cosX;
    result(2, 1) = sinX;
    result(0, 2) = sinY;
    result(1, 2) = -sinX * cosY;
    result(2, 2) = cosX * cosY;

    return result;
  }

  mat4f perspective(float fov, float aspect, float zNear, float zFar)
  {
    assert(std::abs(aspect - std::numeric_limits<float>::epsilon()) > 0.0);

    const float tanHalfFov = std::tan(fov * 0.5);
    mat4f result = mat4f::Zero();
    result(0, 0) = 1.0 / (aspect * tanHalfFov);
    result(1, 1) = 1.0 / (tanHalfFov);
    result(2, 2) = -(zFar + zNear) / (zFar - zNear);
    result(2, 3) = -(2.0 * zFar * zNear) / (zFar - zNear);
    result(3, 2) = -1.0;

    return result;
  }

  mat4f ortho(float left, float right, float bottom, float top, float zNear, float zFar)
  {
    mat4f result = mat4f::Identity();
    result(0, 0) = 2.0 / (right - left);
    result(1, 1) = 2.0 / (top - bottom);
    result(2, 2) = -2.0 / (zFar - zNear);
    result(0, 3) = -(right + left) / (right - left);
    result(1, 3) = -(top + bottom) / (top - bottom);
    result(2, 3) = -(zFar + zNear) / (zFar - zNear);

    return result;
  }

  vec3f center(vec3f const &a, vec3f const &b)
  {
    vec3f ret;
    ret.x() = (a.x() + b.x()) * 0.5f;
    ret.y() = (a.y() + b.y()) * 0.5f;
    ret.z() = (a.z() + b.z()) * 0.5f;

    return ret;
  }

  float distance2(vec3f const &a, vec3f const &b)
  {
    float dx = (b.x() - a.x());
    float dy = (b.y() - a.y());
    float dz = (b.z() - a.z());

    return dx * dx + dy * dy + dz * dz;
  }

  float distance(vec3f const &a, vec3f const &b)
  {
    return sqrt(distance2(a, b));
  }

  vec3f subVec(vec3f const &u, vec3f const &v)
  {
    vec3f w;
    w.x() = v.x() - u.x();
    w.y() = v.y() - u.y();
    w.z() = v.z() - u.z();

    return w;
  }

  vec3f mult_vec_mat(vec3f const &u, mat4f const &m)
  {
    vec4f v = u.homogeneous();
    float x = m.row(0).dot(v);
    float y = m.row(1).dot(v);
    float z = m.row(2).dot(v);
    float w = m.row(3).dot(v);

    assert(w != 0);
    float d = 1.f / w;
    if (w == 1.f)
      return vec3f(x, y, z);
    return vec3f(x * d, y * d, z * d);

    return u;
  }

  vec3f lerp(const vec3f& u, const vec3f& v, const float t)
  {
    return u * (1 - t) + v * t;
  }
} // namespace utils


namespace transform
{
  mat4f rotationX(float theta)
  {
    Eigen::Affine3f transform{Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitX()).toRotationMatrix()};
    return transform.matrix();
  }

  mat4f rotationY(float theta)
  {
    Eigen::Affine3f transform{Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitY()).toRotationMatrix()};
    return transform.matrix();
  }

  mat4f rotationZ(float theta)
  {
    Eigen::Affine3f transform{Eigen::AngleAxisf(theta, Eigen::Vector3f::UnitZ()).toRotationMatrix()};
    return transform.matrix();
  }

  mat4f rotation(float theta, vec3f axis)
  {
    Eigen::Affine3f transform{Eigen::AngleAxisf(theta, axis).toRotationMatrix()};
    return transform.matrix();
  }

  mat4f rotation(const quatf& quaternion)
  {
    mat3f rotation3x3 = quaternion.toRotationMatrix();
    mat4f rotation = mat4f::Identity();
    rotation.block<3, 3>(0, 0) = rotation3x3;

    return rotation;
  }

  mat4f translation(const vec3f& v)
  {
    Eigen::Affine3f transform{Eigen::Translation3f(v)};
    return transform.matrix();
  }

  mat4f translation(const float x, const float y, const float z)
  {
    vec3f v(x, y, z);
    Eigen::Affine3f transform{Eigen::Translation3f(v)};
    return transform.matrix();
  }

  mat4f viewport(const float width, const float height)
  {
    float w = width * .5f;
    float h = height * .5f;

    mat4f ret;
    ret << w, 0, 0, w,
        0, h, 0, h,
        0, 0, .5f, .5f,
        0, 0, 0, 1.f;

    return ret;
  }
}; // transform

#endif // CGAL_UTILS_EIGEN_H
