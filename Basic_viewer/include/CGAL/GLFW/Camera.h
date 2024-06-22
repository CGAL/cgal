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
    Camera();

    void lookat(vec3f const &center, const float size);
    void lookat(vec3f const &pmin, vec3f const &pmax);

    void update(const float dt);

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

    void set_fov(float d);

    inline void set_orientation(const quatf& orientation) { m_orientation = orientation; }

    inline float get_size() const { return m_size; }

    void reset_position();
    void reset_rotation();
    void reset_all();
    
    inline float get_radius() const { return m_radius; }
    inline vec3f get_center() const { return m_center; }

    inline void mode(MODE mode) { m_mode = mode; }
    inline void toggle_mode() { m_mode = m_mode == ORTHOGRAPHIC ? PERSPECTIVE : ORTHOGRAPHIC; }
    inline bool is_orthographic() const { return !m_orbiter; }

    inline void inc_rspeed() { m_rspeed = std::min(m_rspeed+10.f, 500.f); }
    inline void dec_rspeed() { m_rspeed = std::max(m_rspeed-10.f, 50.f); }

    inline void inc_tspeed() { m_tspeed = std::min(m_tspeed+.1f, 20.f); }
    inline void dec_tspeed() { m_tspeed = std::max(m_tspeed-.1f, 0.5f); }

    inline void inc_rotation_smoothness() { m_rotationSmoothFactor = std::max(m_rotationSmoothFactor-.01f, 0.01f); }
    inline void dec_rotation_smoothness() { m_rotationSmoothFactor = std::min(m_rotationSmoothFactor+.01f, 1.f); }

    inline void inc_translation_smoothness() { m_translationSmoothFactor = std::max(m_translationSmoothFactor-.01f, 0.01f); }
    inline void dec_translation_smoothness() { m_translationSmoothFactor = std::min(m_translationSmoothFactor+.01f, 1.f); }

    void set_zoom_smoothness(const float s);

    inline void toggle_fly() { m_orbiter = !m_orbiter; reset_all(); }

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

    /* SMOOTHNESS PARAMS */

    float m_targetSize;

    float m_pitch; 
    float m_yaw; 
    
    float m_targetPitch; 
    float m_targetYaw; 

    float m_rotationSmoothFactor;
    float m_translationSmoothFactor;
    float m_zoomSmoothFactor;

    float m_targetFov;

    vec3f m_targetPosition;

};

/********************DEFINITIONS********************/

Camera::Camera() : 
    m_mode(PERSPECTIVE), 
    m_center(), 
    m_position(vec3f::Zero()), 
    m_targetPosition(vec3f::Zero()), 
    m_orientation(quatf::Identity()), 
    m_tspeed(1.f), m_rspeed(150.f), 
    m_size(5.f), m_cst_size(5.f), 
    m_radius(5.f), m_fov(45.f),
    m_width(1.f), m_height(1.f),
    m_orbiter(true),
    m_pitch(0.f),
    m_yaw(0.f),
    m_targetPitch(0.f),
    m_targetYaw(0.f),
    m_targetSize(5.f),
    m_targetFov(45.f),
    m_zoomSmoothFactor(0.1f),
    m_rotationSmoothFactor(0.1f),
    m_translationSmoothFactor(0.1f)
{

}

inline 
void Camera::lookat(vec3f const &center, const float size)
{
    m_center = center;   // camera focus
    m_size = size;
    m_targetSize = size;
    m_cst_size = size; // allowing reset
    m_radius = size;

    m_orientation = quatf::Identity();
}

inline 
void Camera::lookat(vec3f const &pmin, vec3f const &pmax)
{
    lookat(center(pmin, pmax), distance(pmin, pmax));
}

inline 
void Camera::update(const float dt) 
{
    if (!m_orbiter) {
        if (m_targetPitch > 90.f) 
        {
            m_targetPitch = 90.f;
        }
        else if (m_targetPitch < -90.f) 
        {
            m_targetPitch = -90.f;
        }
    }

    float smoothPitch = m_pitch + m_rotationSmoothFactor * (m_targetPitch - m_pitch);   
    float smoothYaw = m_yaw + m_rotationSmoothFactor * (m_targetYaw - m_yaw);   

    float pitchDelta = radians((smoothPitch - m_pitch) * m_rspeed) * dt;
    float yawDelta = radians((smoothYaw - m_yaw) * m_rspeed) * dt;

    if (m_orbiter) 
    {

        quatf rotX(Eigen::AngleAxisf(pitchDelta, vec3f::UnitX()));
        quatf rotY(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
        m_orientation = rotX * rotY * m_orientation;
    }
    else // free fly
    {
        quatf rotX(Eigen::AngleAxisf(pitchDelta, m_orientation.inverse() * vec3f::UnitX()));
        quatf rotY(Eigen::AngleAxisf(yawDelta, vec3f::UnitY()));
        m_orientation = m_orientation * rotX * rotY;
    }
    m_pitch = smoothPitch;
    m_yaw = smoothYaw;

    m_position += m_translationSmoothFactor * (m_targetPosition - m_position);

    m_size += m_zoomSmoothFactor * (m_targetSize - m_size);

    m_fov += .5f * (m_targetFov - m_fov);
}

inline 
void Camera::translation(const float x, const float y)
{
    float xspeed = x * m_tspeed;
    float yspeed = y * m_tspeed;

    if (m_orbiter) 
    {
        m_targetPosition.x() += m_size * xspeed;
        m_targetPosition.y() += m_size * yspeed;
    } 
    else // free fly  
    {
        vec3f right = (m_orientation.inverse() * vec3f::UnitX()).normalized();
        vec3f up = (m_orientation.inverse() * vec3f::UnitY()).normalized();
        m_targetPosition -= right * m_size * xspeed; 
        m_targetPosition -= up * m_size * yspeed; 
    }
}

inline 
void Camera::rotation(const float x, const float y)
{
    m_targetPitch += y;
    m_targetYaw += x;
}

inline 
void Camera::move(const float z)
{
    if (m_orbiter) 
    {
        m_targetSize -= m_targetSize * z;
        if (m_targetSize < 0.001f)
            m_targetSize = 0.001f;
    }
    else // free fly
    {
        vec3f forward = (m_orientation.inverse() * vec3f::UnitZ()).normalized();
        m_targetPosition -= forward * m_size * z * m_tspeed;
    }
}

inline 
void Camera::move_left(const float dt)
{
    move_right(-dt);
}

inline 
void Camera::move_right(const float dt)
{
    if (m_orbiter) 
    {
        m_targetPosition.x() += m_size * dt * m_tspeed;
    }
    else // free fly
    {
        vec3f right = (m_orientation.inverse() * -vec3f::UnitX()).normalized();
        m_targetPosition -= right * m_size * dt * m_tspeed;
    }
}

inline 
void Camera::move_up(const float dt)
{
    if (m_orbiter) 
    {
        m_targetPosition.y() += m_size * dt * m_tspeed;
    }
    else // free fly
    {
        vec3f up = (m_orientation.inverse() * -vec3f::UnitY()).normalized();
        m_targetPosition -= up * m_size * dt * m_tspeed;
    }
}

inline 
void Camera::move_down(const float dt)
{
    move_up(-dt);
}

inline 
mat4f Camera::view() const
{
    mat3f rotation3x3 = m_orientation.toRotationMatrix();
    mat4f rotation = mat4f::Identity();
    rotation.block<3, 3>(0, 0) = rotation3x3;

    if (m_orbiter) 
    {

        mat4f camera = transform::translation(-m_position.x(), -m_position.y(), -m_size); 
        mat4f model = transform::translation(-m_center.x(), -m_center.y(), -m_center.z()); 

        return  camera * rotation * model;
    } 
    else // free fly  
    {
        mat4f translation = transform::translation(-m_position.x(), -m_position.y(), -m_position.z()-m_size); // translate camera to (0,0,0)

        return rotation * translation;
    }
}

inline 
void Camera::set_fov(const float d) 
{ 

    m_targetFov += d * 2.f; 
    if (m_targetFov > 160.f) m_targetFov = 160.;
    if (m_targetFov < 15.f) m_targetFov = 15.f;

    if (m_mode == PERSPECTIVE)
    {
        float aspectRatio = m_width / m_height;
        float tanHalfFov = tan(radians(m_targetFov) * .5f);
        m_targetSize = m_radius / tanHalfFov;
    }
}

inline 
void Camera::set_zoom_smoothness(const float s)
{
    m_zoomSmoothFactor += s * 0.01; 
    if (m_zoomSmoothFactor > 1.f) m_zoomSmoothFactor = 1.;
    if (m_zoomSmoothFactor < 0.01f) m_zoomSmoothFactor = .01;
}

inline 
mat4f Camera::projection(const float width, const float height)
{
    m_width = width;
    m_height = height;

    return projection();
}

inline 
mat4f Camera::projection() const
{ 
    if (m_mode == ORTHOGRAPHIC) {
        float aspect = m_width / m_height;
        float halfWidth = m_size * aspect * 0.5f;
        float halfHeight = m_size * 0.5f;
        return ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, znear(), zfar());
    }

    return perspective(radians(m_fov), m_width / m_height, znear(), zfar());
}

inline 
float Camera::znear() const
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

inline 
float Camera::zfar() const
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

inline 
mat4f Camera::viewport() const
{
    return transform::viewport(m_width, m_height);
}

inline 
vec3f Camera::position() const
{
    mat4f v = view();
    mat4f vi = v.inverse();

    vec3f o(0.f, 0.f, 0.f);
    return multVecMat(o, vi);
}

inline 
vec3f Camera::forward() const
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
        fw = m_orientation.inverse() * -vec3f::UnitZ();
    }

    return fw.normalized();
}

inline 
void Camera::reset_rotation() 
{ 
    m_pitch = 0.f; 
    m_yaw = 0; 
    m_targetPitch = 0.f; 
    m_targetYaw = 0; 
    m_orientation = quatf::Identity(); 
}


inline void Camera::reset_position()
{
    m_position = vec3f::Zero();
    m_targetPosition = vec3f::Zero(); 
}


inline 
void Camera::reset_all()
{
    reset_position();
    reset_rotation();
    m_size = m_cst_size;
    m_targetSize = m_cst_size;
    m_fov = 45.f;
    m_targetFov = 45.f;
    m_zoomSmoothFactor = .1f;
    m_rotationSmoothFactor = .1f;
    m_translationSmoothFactor = .1f;
}

inline 
mat4f Camera::image_toWorld() const 
{
    mat4f v = view();
    mat4f p = projection();
    mat4f vp = viewport();
    mat4f w2i = vp * p * v;    
    return w2i.inverse();
}


inline 
void Camera::frame(const float z, vec3f &dO, vec3f &dx, vec3f &dy) const
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
