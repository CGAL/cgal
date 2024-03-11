/**
 * @file   ui/gl/MouseAdapter.h
 * @author Gernot Walzl
 * @date   2011-12-22
 */

#ifndef UI_GL_MOUSEADAPTER_H
#define UI_GL_MOUSEADAPTER_H

#include "ui/gl/ptrs.h"
#include "ui/gl/typedefs.h"
#include "util/ptrs.h"
#include <string>

namespace ui { namespace gl {

class MouseAdapter {
public:
    static MouseAdapterSPtr create(CameraSPtr camera);
    virtual ~MouseAdapter();

    /**
     * It will return -1 when the given string is invalid.
     */
    int toButton(const std::string& button);
    void loadConfig(util::ConfigurationSPtr config);

    void pressed(int button, int x, int y);
    void released(int button, int x, int y);
    void dragged(int x, int y);

protected:
    MouseAdapter(CameraSPtr camera);

    CameraSPtr camera_;

    // buttons
    int b_look_;
    int b_rotate_;
    int b_pan_;
    int b_forward_;
    int b_backpedal_;

    // state
    bool mouse_look_;
    bool mouse_rotate_;
    bool mouse_pan_;

    // sensitivity
    double sensitivity_look_;
    double sensitivity_rotate_;
    double sensitivity_pan_;
    double sensitivity_move_;

    vec2i prev_pos_;
};

} }

#endif /* UI_GL_MOUSEADAPTER_H */

