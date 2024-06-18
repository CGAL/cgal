#pragma once

#include <glad/glad.h>
#include <map>
#include <string>

class Shader
{
public:
    Shader() : program(0) {}
    Shader(int program) : program(program) {}

    void destroy()
    {
        glDeleteProgram(program);
    }

    void use()
    {
        glUseProgram(program);
    }

    int getUniform(const std::string &name)
    {
        int loc = m_uniforms[name];
        if (loc != 0)
        {
            return loc;
        }

        loc = glGetUniformLocation(program, name.c_str());
        m_uniforms[name] = loc;
        return loc;
    }

    void setMatrix4f(const std::string &name, GLfloat *data, GLboolean transpose = false)
    {
        glUniformMatrix4fv(getUniform(name), 1, transpose, data);
    }

    void setVec4f(const std::string &name, GLfloat *data)
    {
        glUniform4fv(getUniform(name), 1, data);
    }

    void setFloat(const std::string &name, float data)
    {
        glUniform1f(getUniform(name), data);
    }

    static Shader loadShader(std::string src_vertex, std::string src_fragment, std::string name = "")
    {
        unsigned int vshader = glCreateShader(GL_VERTEX_SHADER);
        const char *source_ = src_vertex.c_str();
        glShaderSource(vshader, 1, &source_, NULL);
        glCompileShader(vshader);
        Shader::checkCompileErrors(vshader, "VERTEX", name);

        source_ = src_fragment.c_str();
        unsigned int fshader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fshader, 1, &source_, NULL);
        glCompileShader(fshader);
        Shader::checkCompileErrors(fshader, "FRAGMENT", name);

        unsigned int program = glCreateProgram();
        glAttachShader(program, vshader);
        glAttachShader(program, fshader);

        glLinkProgram(program);
        Shader::checkCompileErrors(program, "PROGRAM", name);

        glDeleteShader(vshader);
        glDeleteShader(fshader);

        return Shader(program);
    }

private:
    std::unordered_map<std::string, int> m_uniforms;
    int program;

    static void checkCompileErrors(GLuint shader, std::string type, std::string name)
    {
        GLint success;
        GLchar infoLog[1024];

        if (type != "PROGRAM")
        {
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

            if (!success)
            {
                glGetShaderInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "_" << name << "\n"
                          << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
        else
        {
            glGetProgramiv(shader, GL_LINK_STATUS, &success);

            if (!success)
            {
                glGetProgramInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "_" << name << "\n"
                          << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
    }
};