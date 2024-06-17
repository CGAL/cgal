#pragma once

#include "math.h"

class Camera {
public:
    Camera() : m_center(),  m_position(), m_orientation(), m_size(5.f), 
    m_cst_size(5.f), m_radius(5.f), m_width(1.f), m_height(1.f), m_fov(45.f) {}

    void lookat(vec3f const& center, const float size) 
    {
        m_center = center;
        m_position = {0, 0};
        m_orientation = {0, 0};
        m_size = size;
        m_cst_size = size;
        m_radius = size;
    }

    void lookat(vec3f const& pmin, vec3f const& pmax) 
    {
        lookat(center(pmin, pmax), distance(pmin, pmax));
    }

    void translation(const float x, const float y) 
    {
        m_position.x() = m_position.x() - m_size * x;
        m_position.y() = m_position.y() + m_size * y;
    }

    void rotation(const float x, const float y) 
    {
        m_orientation.x() = m_orientation.x() + y;
        m_orientation.y() = m_orientation.y() + x;
    }

    void move(const float z) 
    {
        m_size = m_size - m_size * 0.01F * z;
        // std::cout << "m_size: " << m_size << std::endl;
        if ( m_size < 0.001F ) 
            m_size = 0.001F;
    }

    mat4f view() const 
    {
        mat4f c = transform::translation(-m_position.x(), -m_position.y(), -m_size); // position caméra
        mat4f m = transform::translation(-m_center.x(), -m_center.y(), -m_center.z()); // centre de l'objet observé
        mat4f r = transform::rotationX(m_orientation.x()) * transform::rotationY(m_orientation.y()); 

        return c * r * m; // on position 'c' et oriente 'r' la caméra en fonction de l'objet 'm' observé 
    }

    mat4f projection(const float width, const float height, const float fov)  
    {
        m_fov = fov;
        m_width = width;
        m_height = height;
        
        return projection();
    }

    float znear() const 
    {
        float d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
        return std::max(0.1F, d - 2*m_radius);
    }

    float zfar() const 
    {
        float d = distance(m_center, vec3f(m_position.x(), m_position.y(), m_size));
        // std::cout << "distance (zfar): " << d << std::endl;
        return std::max(100.F, d + 2*m_radius);
    }

    mat4f projection() const 
    {
        return perspective(m_fov, m_width / m_height, znear(), zfar());
    }

    mat4f viewport() const 
    {
        return transform::viewport(m_width, m_height);
    }

    // void frame(const float z, vec3f& dO, vec3f& dx, vec3f& dy) const 
    // {
    //     mat4f v = view();
    //     mat4f p = projection();
    //     mat4f vp = viewport();
    //     mat4f w2i = vp * p * v; // passage du monde vers l'image 
    //     mat4f i2w = w2i.inverse(); // passage de l'image vers le monde 

    //     dO = i2w(vec3f(0, 0, z));
    //     vec3f d1 = i2w(vec3f(1, 0, z)); // récupère la colonne x
    //     vec3f d2 = i2w(vec3f(0, 1, z)); // récupère la colonne y

    //     dx = vec3f(dO, d1); 
    //     dy = vec3f(dO, d2);
    // }

    vec3f position() 
    {
        mat4f v = view();
        mat4f vi = v.inverse();

        vec3f o(0.f, 0.f, 0.f);
        return transform::multMatVec(o, vi);
    }

    void reset_all() 
    {
        reset_position();
        reset_rotation();
        m_size = m_cst_size;
        m_radius = m_cst_size;
    }

    inline void reset_position() { m_position = { 0, 0 }; }
    inline void reset_rotation() { m_orientation = { 0, 0 }; }

    inline float radius( ) const { return m_radius; } 

    void update_center(vec3f const& v) 
    { 
        // vec3f pos1 = position();
        mat4f t = transform::translation(v.x(), v.y(), v.z());
        m_center = transform::multMatVec(m_center, t); 
        // vec3f pos = { 
        //     position().x - pos1.x, 
        //     position().y - pos1.y, 
        //     position().z - pos1.z   
        // };

        // m_position = { pos.x, pos.y };


        // std::cout << "cam position : " 
        //         << position().x << " " 
        //         << position().y << " " 
        //         << position().z << std::endl;
        // m_position = { position().x, position().y };
    }

private:
    vec3f m_center; // centre de l'objet observé
    vec2f m_position; // position de la caméra (translations appliquées)
    vec2f m_orientation; // orientation de la caméra (rotations appliquées) 
    float m_size;
    float m_cst_size;
    float m_radius;

    float m_width; // largeur de l'image
    float m_height; // hauteur de l'image
    float m_fov; // FOV pour la projection 
};
