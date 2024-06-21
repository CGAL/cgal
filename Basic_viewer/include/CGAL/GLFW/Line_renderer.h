#ifndef LINE_RENDERER_H
#define LINE_RENDERER_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <vector>

#include "Shader.h"
#include "Line.h"

class Line_renderer {
public:
    Line_renderer() : m_shader(nullptr), m_mvp(mat4f::Identity()), m_index_count(0) {}

    ~Line_renderer() 
    {
        if (m_vao != 0) glDeleteVertexArrays(1, &m_vao);
        if (m_vbo != 0) glDeleteBuffers(1, &m_vbo);
        if (m_ebo != 0) glDeleteBuffers(1, &m_ebo);
    }

    inline void set_mvp(const mat4f& mvp) { m_mvp = mvp; }

    inline void attach_shader(Shader* shader) 
    { 
        m_shader = shader; 
        
        glGenVertexArrays(1, &m_vao);
        glGenBuffers(1, &m_vbo);
        glGenBuffers(1, &m_ebo);

        glBindVertexArray(m_vao);
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);    
    }

    inline void add_line(const Line& line) 
    {
        m_lines.push_back(line);
        const auto& data = line.get_data();
        m_vertices.insert(m_vertices.end(), data.begin(), data.end());
        m_indices.push_back(m_index_count);
        m_indices.push_back(m_index_count+1);
        m_index_count += 2;
    }

    inline void draw()  
    {
        assert(m_shader != nullptr);

        glBindVertexArray(m_vao);

        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * m_vertices.size(), m_vertices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * m_indices.size(), m_indices.data(), GL_STATIC_DRAW);

        m_shader->use();
        for (int i = 0; i < m_lines.size(); ++i) {
            const auto& line = m_lines[i];
            m_shader->setMatrix4f("mvp_matrix", m_mvp.data()); 
            m_shader->setVec3f("color", line.get_color().data());

            glLineWidth(line.get_width());
            glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, (void*)(i * 2 * sizeof(unsigned int)));
        }

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glLineWidth(1.f);
    }

private: 
    Shader* m_shader;
    
    std::vector<Line> m_lines;
    std::vector<float> m_vertices;
    std::vector<unsigned int> m_indices;

    mat4f m_mvp;

    unsigned int m_vao, m_vbo, m_ebo;
    unsigned int m_index_count;
};

#endif // LINE_RENDERER_H



