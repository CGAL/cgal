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
#include <QDebug>

namespace GlSplat {

void Shader::define(const char* name, const char* value)
{
    mDefines[std::string(name)] = value;
}
//--------------------------------------------------------------------------------
bool Shader::loadSources(const char* vsrc, const char* fsrc)
{
    bool allIsOk = true;
    mProgramID = viewer->glCreateProgram();

    std::string defineStr = "#extension GL_ARB_texture_rectangle : enable\n";
    for(DefineMap::iterator it = mDefines.begin() ; it!=mDefines.end() ; ++it)
    {
      defineStr += "#define " + it->first + " " + it->second + "\n";
    }

    if(vsrc)
    {
        GLuint shaderID = viewer->glCreateShader(GL_VERTEX_SHADER);

        std::string source = defineStr + std::string(vsrc);
        const GLchar * arbSource = source.c_str();

        viewer->glShaderSource(shaderID, 1, (const GLchar **)&arbSource, 0);
        viewer->glCompileShader(shaderID);

        int compiled;
        viewer->glGetShaderiv(shaderID,GL_COMPILE_STATUS,&compiled);
        allIsOk = allIsOk && compiled;
        //printInfoLog(shaderID);

        viewer->glAttachShader(mProgramID, shaderID);
    }

    if(fsrc)
    {
        GLuint shaderID = viewer->glCreateShader(GL_FRAGMENT_SHADER);

        std::string source = defineStr + std::string(fsrc);
        const GLchar * arbSource = source.c_str();

        viewer->glShaderSource(shaderID, 1, (const GLchar **)&arbSource, 0);
        viewer->glCompileShader(shaderID);

        int compiled;
        viewer->glGetShaderiv(shaderID,GL_COMPILE_STATUS,&compiled);
        allIsOk = allIsOk && compiled;
        //printInfoLog(shaderID);

        viewer->glAttachShader(mProgramID, shaderID);
    }

    viewer->glLinkProgram(mProgramID);

    int isLinked;
    viewer->glGetProgramiv(mProgramID, GL_LINK_STATUS, &isLinked);
    allIsOk = allIsOk && isLinked;
    mIsValid = isLinked == GL_TRUE;
    printInfoLog(mProgramID);

    return allIsOk;
}
//--------------------------------------------------------------------------------
void Shader::activate()
{
    assert(mIsValid);
    viewer->glUseProgram(mProgramID);
}
void Shader::release(void)
{
    viewer->glUseProgram(0);
}
//--------------------------------------------------------------------------------
int Shader::getUniformLocation(const char* name)
{
    assert(mIsValid);
    int loc = viewer->glGetUniformLocation(mProgramID, name);
    return loc;
}
//--------------------------------------------------------------------------------
void Shader::setSamplerUnit(const char* sampler, int unit)
{
    activate();
    viewer->glUniform1i(getUniformLocation(sampler), unit);
    release();
}
//--------------------------------------------------------------------------------
int Shader::getAttribLocation(const char* name)
{
    assert(mIsValid);
    int loc = viewer->glGetAttribLocation(mProgramID, name);
    return loc;
}
//--------------------------------------------------------------------------------
void Shader::printInfoLog(GLuint objectID)
{
    int infologLength, charsWritten;
    GLchar *infoLog;
    viewer->glGetProgramiv(objectID,GL_INFO_LOG_LENGTH, &infologLength);
    if(infologLength > 0)
    {
        infoLog = new GLchar[infologLength];
        viewer->glGetProgramInfoLog(objectID, infologLength, &charsWritten, infoLog);
        if (charsWritten>0)
          std::cerr << "Shader info : \n" << infoLog << std::endl;
        delete[] infoLog;
    }
}

void Shader::setViewer(Viewer_interface *v)
{
    viewer = v;
}
} // namepsace GlSplat
