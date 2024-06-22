#pragma once

#include "math.h"

class Camera
{
public:
    enum MODE
    {
        PERSPECTIVE,
        ORTHOGRAPHIC,
        NB_ELT
    };

public:
    Camera() 
    : m_mode(PERSPECTIVE), 
      m_center(), 
      m_position(vec3f::Zero()), 
      m_orientation(quatf::Identity()), 
      m_tspeed(1.f), m_rspeed(5.f), 
      m_size(5.f), m_cst_size(5.f), 
      m_radius(5.f), m_fov(45.f),
      m_width(1.f), m_height(1.f),
      m_orbiter(true) {}

    void lookat(vec3f const &center, const float size);
    void lookat(vec3f const &pmin, vec3f const &pmax);

    void translation(const float x, const float y);
    void rotation(const float x, const float y);
    void move(const float z);

    void move_up(const float dt);
    void move_down(const float dt);
    void move_right(const float dt);
    void move_left(const float dt);

    mat4f view() const;

    mat4f projection() const;
    mat4f projection(const float width, const float height);

    float znear() const;
    float zfar() const;

    mat4f viewport() const;

    mat4f image_toWorld() const;
    void frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const;

    vec3f position() const;
    vec3f forward() const;

    void reset_all();

    void modify_fov(float d);

    inline void set_orientation(const quatf& orientation) { m_orientation = orientation; }

    inline float get_size() const { return m_size; }

    inline void reset_position() { m_position = {0, 0, 0}; }
    inline void reset_rotation() { m_orientation = quatf::Identity(); }

    inline float get_radius() const { return m_radius; }
    inline vec3f get_center() const { return m_center; }

    inline void mode(MODE mode) { m_mode = mode; }
    inline void toggle_mode() { m_mode = m_mode == ORTHOGRAPHIC ? PERSPECTIVE : ORTHOGRAPHIC; }
    inline bool is_orthographic() const { return !m_orbiter; }

    inline void inc_rspeed() { m_rspeed = std::min(m_rspeed+1.f, 20.f); std::cout << m_rspeed << std::endl; }
    inline void dec_rspeed() { m_rspeed = std::max(m_rspeed-1.f, 1.f); std::cout << m_rspeed << std::endl; }

    inline void inc_tspeed() { m_tspeed = std::min(m_tspeed+.1f, 20.f); }
    inline void dec_tspeed() { m_tspeed = std::max(m_tspeed-.1f, 0.5f); }

    inline void toggle_fly() { m_orbiter = !m_orbiter; }



private:
    vec3f m_center;   // centre de l'objet observé
    vec3f m_position; // position de la caméra (translations appliquées)

    quatf m_orientation; // orientation de la caméra (rotations appliquées)

    float m_size;
    float m_cst_size;
    float m_radius;

    float m_width;  // largeur de l'image
    float m_height; // hauteur de l'image
    float m_fov;    // FOV pour la projection
    
    float m_tspeed; 
    float m_rspeed; 

    bool m_orbiter;

    MODE m_mode;
};

/********************DEFINITIONS********************/

inline void Camera::lookat(vec3f const &center, const float size)
{
    m_center = center;   // camera focus
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
    if (m_orbiter) 
    {
        m_position.x() = m_position.x() - m_size * x * m_tspeed;
        m_position.y() = m_position.y() - m_size * y * m_tspeed;
    } 
    else // free fly  
    {
        vec3f right = vec3f::UnitX();
        vec3f up = vec3f::UnitY();
        m_position += right * x * m_tspeed; 
        m_position += up * y * m_tspeed; 
    }
}

inline void Camera::rotation(const float x, const float y)
{
    float pitchDelta = y * m_rspeed;
    float yawDelta = x * m_rspeed;
    quatf rotX(Eigen::AngleAxisf(pitchDelta, vec3f::UnitX()));
    quatf rotY(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
    if (m_orbiter) 
    {
        m_orientation = rotX * rotY * m_orientation;
    }
    else // free fly
    {
        m_orientation = rotY * m_orientation * rotX;
    }
}

inline void Camera::move(const float z)
{
    if (m_orbiter) 
    {
        m_size = m_size - m_size * z;
        // std::cout << "m_size: " << m_size << std::endl;
        if (m_size < 0.001f)
            m_size = 0.001f;
    }
    else // free fly
    {
        vec3f forward = (m_orientation * -vec3f::UnitZ()).normalized();
        vec3f right = (m_orientation * vec3f::UnitX()).normalized();
        vec3f up = (m_orientation * vec3f::UnitY()).normalized();

        std::cout << "orientation : " 
            << m_orientation.x() << " " 
            << m_orientation.y() << " " 
            << m_orientation.z() << " " 
            << m_orientation.w() << " \nfw : (" 
            << forward.x() << ", " 
            << forward.y() << ", " 
            << forward.z() << ") \nrt : ("
            << right.x() << ", " 
            << right.y() << ", " 
            << right.z() << ") \nup : ("
            << up.x() << ", " 
            << up.y() << ", " 
            << up.z() << ")" <<      
        std::endl; 
        m_position += forward * m_size * z * m_tspeed;
    }
}

inline void Camera::move_left(const float dt)
{
    move_right(-dt);
}

inline void Camera::move_right(const float dt)
{
    if (m_orbiter) 
    {
        m_position.x() += m_size * dt * m_tspeed;
    }
    else // free fly
    {
        vec3f right = vec3f::UnitX();
        m_position += right * m_size * dt * m_tspeed;
    }
}

inline void Camera::move_up(const float dt)
{
    if (m_orbiter) 
    {
        m_position.y() += m_size * dt * m_tspeed;
    }
    else // free fly
    {
        vec3f up = vec3f::UnitY();
        m_position += up * m_size * dt * m_tspeed;
    }
}

inline void Camera::move_down(const float dt)
{
    move_up(-dt);
}

inline mat4f Camera::view() const
{
    mat3f rotation3x3 = m_orientation.toRotationMatrix();
    mat4f rotation = mat4f::Identity();
    rotation.block<3, 3>(0, 0) = rotation3x3;

    if (m_orbiter) 
    {

        mat4f position = transform::translation(m_position.x(), m_position.y(), -m_size); // translate camera to (0,0,0)
        mat4f model = transform::translation(-m_center.x(), -m_center.y(), -m_center.z());  // translate focus to (0,0,0)

        return  position * rotation * model;
    } 
    else // free fly  
    {
        mat4f position = transform::translation(m_position.x(), m_position.y(), m_position.z()-m_size); // translate camera to (0,0,0)

        return rotation * position;
    }
}

inline void Camera::modify_fov(float d) { 
    m_fov += d; 
    if (m_fov > 179.f) m_fov = 179.f;
    if (m_fov < 1.f) m_fov = 1.f;
    std::cout << m_fov << std::endl;
}

inline mat4f Camera::projection(const float width, const float height)
{
    m_width = width;
    m_height = height;

    return projection();
}

inline mat4f Camera::projection() const
{
    if (m_mode == ORTHOGRAPHIC) {
        float aspect = m_width / m_height;
        float halfWidth = m_size * aspect * 0.5f;
        float halfHeight = m_size * 0.5f;
        return ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, znear(), zfar());
    }

    return perspective(m_fov, m_width / m_height, znear(), zfar());
}

inline float Camera::znear() const
{
    float d;
    if (m_orbiter) 
    {
        d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
    }
    else
    {
        d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
    }
    return std::max(0.1f, d - 2 * m_radius);
}

inline float Camera::zfar() const
{
    float d;
    if (m_orbiter) 
    {
        d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
    }
    else
    {
        d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_position.z()+m_size));
    }
    return std::max(100.f, d + 2 * m_radius);
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
    if (m_orbiter) 
    {
        fw.x() = m_center.x() - position().x();
        fw.y() = m_center.y() - position().y();
        fw.z() = m_center.z() - position().z();
    }
    else 
    {
        fw = m_orientation * -vec3f::UnitZ();
    }

    return fw.normalized();
}

inline void Camera::reset_all()
{
    reset_position();
    reset_rotation();
    m_size = m_cst_size;
    m_radius = m_cst_size;
}

inline mat4f Camera::image_toWorld() const 
{
    mat4f v = view();
    mat4f p = projection();
    mat4f vp = viewport();
    mat4f w2i = vp * p * v;    
    return w2i.inverse();
}


inline void Camera::frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const
{
    mat4f i2w = image_toWorld(); // passage de l'image vers le monde

    vec3f uo(0, 0, z);
    vec3f ux(1, 0, z);
    vec3f uy(0, 1, z);
    dO = multVecMat(uo, i2w);
    vec3f d1 = multVecMat(ux, i2w); // récupère la colonne x
    vec3f d2 = multVecMat(uy, i2w); // récupère la colonne y

    dx = subVec(dO, d1);
    dy = subVec(dO, d2);
}
