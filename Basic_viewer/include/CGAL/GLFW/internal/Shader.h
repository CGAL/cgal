#ifndef CGAL_SHADER_GLFW_H
#define CGAL_SHADER_GLFW_H

#include <glad/glad.h>
#include <map>
#include <string>

class Shader
{
public:
  Shader() : program(0) {}
  Shader(int program) : program(program) {}

  inline 
  void destroy()
  {
    glDeleteProgram(program);
  }

  inline 
  void use()
  {
    glUseProgram(program);
  }

  inline 
  int get_uniform(const std::string &name)
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

  inline 
  void set_mat4f(const std::string &name, const GLfloat* data, GLboolean transpose = false)
  {
    glUniformMatrix4fv(get_uniform(name), 1, transpose, data);
  }

  inline 
  void set_vec4f(const std::string &name, const GLfloat* data)
  {
    glUniform4fv(get_uniform(name), 1, data);
  }

  inline 
  void set_vec3f(const std::string &name, const GLfloat* data)
  {
    glUniform3fv(get_uniform(name), 1, data);
  }

  inline 
  void set_float(const std::string &name, const float data)
  {
    glUniform1f(get_uniform(name), data);
  }

  inline 
  static Shader load_shader(std::string srcVertex, std::string srcFragment, std::string name = "")
  {
    unsigned int vshader = glCreateShader(GL_VERTEX_SHADER);
    const char* SOURCE = srcVertex.c_str();
    glShaderSource(vshader, 1, &SOURCE, NULL);
    glCompileShader(vshader);
    Shader::check_compile_errors(vshader, "VERTEX", name);

    SOURCE = srcFragment.c_str();
    unsigned int fshader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fshader, 1, &SOURCE, NULL);
    glCompileShader(fshader);
    Shader::check_compile_errors(fshader, "FRAGMENT", name);

    unsigned int program = glCreateProgram();
    glAttachShader(program, vshader);
    glAttachShader(program, fshader);

    glLinkProgram(program);
    Shader::check_compile_errors(program, "PROGRAM", name);

    glDeleteShader(vshader);
    glDeleteShader(fshader);

    return Shader(program);
  }

private:
  std::unordered_map<std::string, int> m_uniforms;
  int program;

  inline 
  static void check_compile_errors(GLuint shader, std::string type, std::string name)
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


#endif // CGAL_SHADER_GLFW_H