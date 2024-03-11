/**
 * @file   ui/gl/Camera.h
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#ifndef UI_GL_CAMERA_H
#define UI_GL_CAMERA_H

#include "ui/gl/ptrs.h"
#include "util/ptrs.h"
#include <GL/gl.h>
#include <string>

namespace ui { namespace gl {

/*!
 * Quake style camera movement
 */
class Camera {

public:
    virtual ~Camera();
    static CameraSPtr create();

    void loadConfig(util::ConfigurationSPtr config);

    void reset();
    void topdown();

    virtual std::string toString();

    GLdouble eye(unsigned int i);
    GLdouble center(unsigned int i);
    GLdouble up(unsigned int i);

    void move(double dx, double dy, double dz);
    void moveFB(double amount);    // forward / backpedal
    void moveUD(double amount);    // up / down  (z-direction)
    void strafeLR(double amount);  // left / right
    void strafeUD(double amount);  // up / down  (orientation dependent)

    void lookLR(double amount);  // left / right
    void lookUD(double amount);  // up / down

    void rotateWorldLR(double amount);  // left / right
    void rotateWorldUD(double amount);  // up / down

protected:
    Camera();

    GLdouble eye_[3];     // current position
    GLdouble center_[3];  // look at
    GLdouble up_[3];      // up vector

    double angleXY();
    double angleZ();
};

} }

#endif /* UI_GL_CAMERA_H */
