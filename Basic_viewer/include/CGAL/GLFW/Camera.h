#pragma once

#include "math.h"

class Camera
{
public:
    enum MODE
    {
        ORBIT,
        FLY,
        ORTHO,
        NB_ELT
    };

public:
    Camera() : m_mode(ORBIT), m_center(), m_position(), m_orientation(quatf::Identity()), m_size(5.f),
               m_cst_size(5.f), m_radius(5.f), m_width(1.f), m_height(1.f), m_fov(45.f) {}

    inline void lookat(vec3f const &center, const float size);
    inline void lookat(vec3f const &pmin, vec3f const &pmax);

    inline void translation(const float x, const float y);
    inline void rotation(const float x, const float y);
    inline void move(const float z);

    inline mat4f view() const;
    inline mat4f projection(const float width, const float height, const float fov);

    inline float znear() const;
    inline float zfar() const;

    inline mat4f projection() const;

    inline mat4f viewport() const;

    inline void frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const;

    inline vec3f position() const;
    inline vec3f forward() const;

    inline void reset_all();

    inline void reset_position() { m_position = {0, 0}; }
    inline void reset_rotation() { m_orientation = quatf::Identity(); }

    inline float get_radius() const { return m_radius; }
    inline vec3f get_center() const { return m_center; }

    inline void mode(MODE mode) { m_mode = mode; }

private:
    inline mat4f view_orbit() const;
    inline mat4f view_fly() const;
    inline mat4f view_ortho() const;
    inline mat4f projection_orbit() const;
    inline mat4f projection_fly() const;
    inline mat4f projection_ortho() const;

private:
    vec3f m_center;   // centre de l'objet observé
    vec2f m_position; // position de la caméra (translations appliquées)

    quatf m_orientation; // orientation de la caméra (rotations appliquées)

    float m_size;
    float m_cst_size;
    float m_radius;

    float m_width;  // largeur de l'image
    float m_height; // hauteur de l'image
    float m_fov;    // FOV pour la projection

    MODE m_mode;
};

/********************DEFINITIONS********************/

inline void Camera::lookat(vec3f const &center, const float size)
{
    m_center = center;   // camera focus
    m_position = {0, 0}; // camera position
    m_size = size;
    m_cst_size = size; // allowing reset
    m_radius = size;

    m_orientation = quatf::Identity();
}

inline void Camera::lookat(vec3f const &pmin, vec3f const &pmax)
{
    lookat(center(pmin, pmax), distance(pmin, pmax));
}

inline void Camera::translation(const float x, const float y)
{
    m_position.x() = m_position.x() - m_size * x;
    m_position.y() = m_position.y() + m_size * y;
}

inline void Camera::rotation(const float x, const float y)
{
    float alpha = radians(x);
    float beta = radians(y);
    quatf rotX(Eigen::AngleAxisf(beta, vec3f::UnitX()));
    quatf rotY(Eigen::AngleAxisf(alpha, vec3f::UnitY()));
    m_orientation = rotY * m_orientation * rotX;
}

inline void Camera::move(const float z)
{
    m_size = m_size - m_size * 0.01f * z;
    // std::cout << "m_size: " << m_size << std::endl;
    if (m_size < 0.001f)
        m_size = 0.001f;
}

inline mat4f Camera::view() const
{
    mat4f position = transform::translation(-m_position.x(), -m_position.y(), -m_size); // translate camera to (0,0,0)
    mat4f model = transform::translation(-m_center.x(), -m_center.y(), -m_center.z());  // translate focus to (0,0,0)

    mat3f rotation3x3 = m_orientation.toRotationMatrix();
    mat4f rotation = mat4f::Identity();
    rotation.block<3, 3>(0, 0) = rotation3x3;

    return position * rotation * model; // on positionne 'position' et oriente 'rotation' la caméra en fonction de l'objet 'model' observé
}

inline mat4f Camera::projection(const float width, const float height, const float fov)
{
    m_fov = fov;
    m_width = width;
    m_height = height;

    return projection();
}

inline float Camera::znear() const
{
    float d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
    return std::max(0.1f, d - 2 * m_radius);
}

inline float Camera::zfar() const
{
    float d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
    // std::cout << "distance (zfar): " << d << std::endl;
    return std::max(100.f, d + 2 * m_radius);
}

inline mat4f Camera::projection() const
{
    return perspective(m_fov, m_width / m_height, znear(), zfar());
}

inline mat4f Camera::viewport() const
{
    return transform::viewport(m_width, m_height);
}

inline vec3f Camera::position() const
{
    mat4f v = view();
    mat4f vi = v.inverse();

    vec3f o(0.f, 0.f, 0.f);
    return multVecMat(o, vi);
}

inline vec3f Camera::forward() const
{
    vec3f fw;
    fw.x() = m_center.x() - position().x();
    fw.y() = m_center.y() - position().y();
    fw.z() = m_center.z() - position().z();

    return fw.normalized();
}

inline void Camera::reset_all()
{
    reset_position();
    reset_rotation();
    m_size = m_cst_size;
    m_radius = m_cst_size;
}

inline void Camera::frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const
{
    mat4f v = view();
    mat4f p = projection();
    mat4f vp = viewport();
    mat4f w2i = vp * p * v;    // passage du monde vers l'image
    mat4f i2w = w2i.inverse(); // passage de l'image vers le monde

    vec3f uo(0, 0, z);
    vec3f ux(1, 0, z);
    vec3f uy(0, 1, z);
    dO = multVecMat(uo, i2w);
    vec3f d1 = multVecMat(ux, i2w); // récupère la colonne x
    vec3f d2 = multVecMat(uy, i2w); // récupère la colonne y

    dx = subVec(dO, d1);
    dy = subVec(dO, d2);
}

inline mat4f Camera::view_orbit() const
{
    return mat4f::Identity();
}

inline mat4f Camera::view_fly() const
{
    return mat4f::Identity();
}

inline mat4f Camera::view_ortho() const
{
    return mat4f::Identity();
}

inline mat4f Camera::projection_orbit() const
{
    return mat4f::Identity();
}

inline mat4f Camera::projection_fly() const
{
    return mat4f::Identity();
}

inline mat4f Camera::projection_ortho() const
{
    return mat4f::Identity();
}
