// This file is part of GlSplat, a simple splatting C++ library
//
// Copyright (C) 2008-2009 Gael Guennebaud <g.gael@free.fr>
//
// GlSplat is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// GlSplat is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with GlSplat. If not, see <http://www.gnu.org/licenses/>.

#include "Shader.h"
#include <iostream>

#include <qfile.h>
#include <qtextstream.h>
#include <assert.h>

namespace GlSplat {

void Shader::define(const char* name, const char* value)
{
    mDefines[std::string(name)] = value;
}
//--------------------------------------------------------------------------------
bool Shader::loadSources(const char* vsrc, const char* fsrc)
{
    bool allIsOk = false;

    mProgramID = glCreateProgram();

    std::string defineStr = "";
    for(DefineMap::iterator it = mDefines.begin() ; it!=mDefines.end() ; ++it)
    {
      defineStr += "#define " + it->first + " " + it->second + "\n";
    }

    if(vsrc)
    {
        GLuint shaderID = glCreateShader(GL_VERTEX_SHADER);

        std::string source = defineStr + std::string(vsrc);
        const GLchar * arbSource = source.c_str();

        glShaderSource(shaderID, 1, (const GLchar **)&arbSource, 0);
        glCompileShader(shaderID);

        int compiled;
        glGetShaderiv(shaderID,GL_COMPILE_STATUS,&compiled);
        allIsOk = allIsOk && compiled;
        //printInfoLog(shaderID);

        glAttachShader(mProgramID, shaderID);
    }

    if(fsrc)
    {
        GLuint shaderID = glCreateShader(GL_FRAGMENT_SHADER);

        std::string source = defineStr + std::string(fsrc);
        const GLchar * arbSource = source.c_str();

        glShaderSource(shaderID, 1, (const GLchar **)&arbSource, 0);
        glCompileShader(shaderID);

        int compiled;
        glGetShaderiv(shaderID,GL_COMPILE_STATUS,&compiled);
        allIsOk = allIsOk && compiled;
        //printInfoLog(shaderID);

        glAttachShader(mProgramID, shaderID);
    }

    glLinkProgram(mProgramID);

    int isLinked;
    glGetProgramiv(mProgramID, GL_LINK_STATUS, &isLinked);
    allIsOk = allIsOk && isLinked;
    mIsValid = isLinked == GL_TRUE;
    printInfoLog(mProgramID);

    return allIsOk;
}
//--------------------------------------------------------------------------------
void Shader::activate(void) const
{
    assert(mIsValid);
    glUseProgram(mProgramID);
}
void Shader::release(void) const
{
    glUseProgram(0);
}
//--------------------------------------------------------------------------------
int Shader::getUniformLocation(const char* name) const
{
    assert(mIsValid);
    int loc = glGetUniformLocation(mProgramID, name);
    return loc;
}
//--------------------------------------------------------------------------------
void Shader::setSamplerUnit(const char* sampler, int unit) const
{
    activate();
    glUniform1i(getUniformLocation(sampler), unit);
    release();
}
//--------------------------------------------------------------------------------
int Shader::getAttribLocation(const char* name) const
{
    assert(mIsValid);
    int loc = glGetAttribLocation(mProgramID, name);
    return loc;
}
//--------------------------------------------------------------------------------
void Shader::printInfoLog(GLuint objectID)
{
    int infologLength, charsWritten;
    GLchar *infoLog;
    glGetProgramiv(objectID,GL_INFO_LOG_LENGTH, &infologLength);
    if(infologLength > 0)
    {
        infoLog = new GLchar[infologLength];
        glGetProgramInfoLog(objectID, infologLength, &charsWritten, infoLog);
        if (charsWritten>0)
          std::cerr << "Shader info : \n" << infoLog << std::endl;
        delete[] infoLog;
    }
}

} // namepsace GlSplat
