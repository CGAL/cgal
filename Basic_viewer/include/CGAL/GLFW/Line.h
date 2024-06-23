#ifndef LINE_H
#define LINE_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>

#include "Shader.h"
#include "math.h"

// https://stackoverflow.com/questions/73588117/opengl-how-to-draw-a-axis-indicator-at-the-top-right-corner (thanks)

class Line {
public:
    Line(const vec3f& start, const vec3f& end, const vec3f& color, const float width) 
    : m_start(start), m_end(end), m_color(color), m_width(width)
    {
        m_data = { 
            start.x(), start.y(), start.z(), 1.f, 
            end.x()  , end.y()  , end.z()  , 1.f 
        };
    }
 
    inline void set_color(vec3f color) { m_color = color; }
    inline void get_width(const float width) { m_width = width; }

    inline const vec3f& get_color() const { return m_color; }
    inline const float& get_width() const { return m_width; }
    inline const std::vector<float>& get_data() const { return m_data; }

private: 
    std::vector<float> m_data;
    vec3f m_start;
    vec3f m_end;
    vec3f m_color;
    float m_width;
};

#endif // LINE_H