/**
 * @file   ui/gl/KeyboardAdapter.h
 * @author Gernot Walzl
 * @date   2011-12-22
 */

#ifndef UI_GL_KEYBOARDADAPTER_H
#define UI_GL_KEYBOARDADAPTER_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "util/ptrs.h"
#include "ui/gl/ptrs.h"
#include <string>

namespace ui { namespace gl {

using algo::ControllerSPtr;

class KeyboardAdapter {
public:
    static KeyboardAdapterSPtr create(CameraSPtr camera, ControllerSPtr controller);
    virtual ~KeyboardAdapter();

    void setWindow(MainOpenGLWindowSPtr window);

    /**
     * special keys (glutSpecialFunc(...)) have negative values
     */
    int toKey(const std::string& key);
    void loadConfig(util::ConfigurationSPtr config);

    void pressed(int key, int x, int y);
    void released(int key, int x, int y);

protected:
    KeyboardAdapter(CameraSPtr camera, ControllerSPtr controller);

    void run();  // fixed rate timer

    ThreadSPtr thread_;
    CameraSPtr camera_;
    ControllerSPtr controller_;
    MainOpenGLWindowWPtr window_;

    // keys
    int k_move_forward_;
    int k_move_backpedal_;
    int k_move_left_;
    int k_move_right_;
    int k_move_up_;
    int k_move_down_;
    int k_look_left_;
    int k_look_right_;
    int k_look_up_;
    int k_look_down_;
    int k_reset_;
    int k_pause_;
    int k_step_;
    int k_skip_;
    int k_toggle_roof_;
    int k_toggle_poly_;
    int k_toggle_skel_;
    int k_toggle_experimental_;
    int k_inc_thickness_;
    int k_dec_thickness_;
    int k_save_skel_;
    int k_save_last_poly_;
    int k_dump_window_;
    int k_print_screen_;
    int k_print_cutpattern_;

    // state
    bool move_forward_;
    bool move_backpedal_;
    bool move_left_;
    bool move_right_;
    bool move_up_;
    bool move_down_;
    bool look_left_;
    bool look_right_;
    bool look_up_;
    bool look_down_;

    double speed_move_;
    double speed_look_;
};

} }

#endif /* UI_GL_KEYBOARDADAPTER_H */
