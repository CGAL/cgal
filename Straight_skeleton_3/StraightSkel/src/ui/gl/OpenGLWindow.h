/**
 * @file   ui/gl/OpenGLWindow.h
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#ifndef UI_GL_OPENGLWINDOW_H
#define UI_GL_OPENGLWINDOW_H

#include "typedefs_thread.h"
#include "ui/gl/ptrs.h"
#include "ui/gl/typedefs.h"
#include <GL/glut.h>
#include <map>
#include <string>

namespace ui { namespace gl {

/*!
 * This class is intended do be subclassed.
 */
class OpenGLWindow : public std::enable_shared_from_this<OpenGLWindow> {

public:
    virtual ~OpenGLWindow();
    static OpenGLWindowSPtr create(int argc, const char* argv[],
            unsigned int width, unsigned int height,
            const char* title);
    void createWindow();
    static void mainLoop();
    bool writeBMP(const char* filename, unsigned int width, unsigned int height, unsigned char* image);
    bool dumpWindow(const char* filename);

protected:
    static std::map<int, OpenGLWindowWPtr> windows_;
    static bool glut_started_;
    static RecursiveMutex mutex_;

    OpenGLWindow(int argc, const char* argv[],
            unsigned int width, unsigned int height,
            const char* title);

    int argc_;
    const char** argv_;
    const char* title_;
    unsigned int width_;
    unsigned int height_;
    double fovy_;
    int window_;
    CameraSPtr camera_;
    vec4f color_;

    void initLighting();
    void initBlending();
    void init();

    // drawing functions
    void setColor(const vec4f rgba);
    void getColor(vec4f& out) const;
    void drawSphere(const vec3f position, float radius);
    void drawCylinder(const vec3f src, const vec3f dst, float radius);
    void drawCircularCylinder(const vec3f center, const vec3f axis,
            const vec3f src, const vec3f dst, float radius);
    void drawTriangle(const vec3f a, const vec3f b, const vec3f c);
    void drawArrow(const vec3f position, const vec3f direction);
    void drawText(const std::string& text);
    void drawCrosshair(float size);

    // to be implemented by inherited class
    virtual void drawContent();
    virtual void handleKeyPressed(int key, int x, int y);
    virtual void handleKeyReleased(int key, int x, int y);
    virtual void handleMousePressed(int button, int x, int y);
    virtual void handleMouseReleased(int button, int x, int y);
    virtual void handleMotion(int x, int y);

    // static functions for glut
    static void displayFunc();
    static void reshapeFunc(int width, int height);
    static void keyFunc(unsigned char key, int x, int y);
    static void keyUpFunc(unsigned char key, int x, int y);
    static void specialFunc(int key, int x, int y);
    static void specialUpFunc(int key, int x, int y);
    static void mouseFunc(int button, int state, int x, int y);
    static void motionFunc(int x, int y);
};

} }

#endif /* UI_GL_OPENGLWINDOW_H */
