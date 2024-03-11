/**
 * @file   ui/gl/KeyboardAdapter.cpp
 * @author Gernot Walzl
 * @date   2011-12-22
 */

#include "ui/gl/KeyboardAdapter.h"

#include "debug.h"
#include "algo/Controller.h"
#include "util/Configuration.h"
#include "ui/gl/Camera.h"
#include "ui/gl/MainOpenGLWindow.h"

namespace ui { namespace gl {

KeyboardAdapter::KeyboardAdapter(CameraSPtr camera, ControllerSPtr controller) {
    this->camera_ = camera;
    this->controller_ = controller;

    this->k_move_forward_ = 'w';
    this->k_move_backpedal_ = 's';
    this->k_move_left_ = 'a';
    this->k_move_right_ = 'd';
    this->k_move_up_ = 'e';
    this->k_move_down_ = 'c';
    this->k_look_left_ = -100;
    this->k_look_right_ = -102;
    this->k_look_up_ = -101;
    this->k_look_down_ = -103;
    this->k_reset_ = 8;
    this->k_pause_ = 'p';
    this->k_step_ = ' ';
    this->k_skip_ = '\033';
    this->k_toggle_roof_ = 'r';
    this->k_toggle_poly_ = 't';
    this->k_toggle_skel_ = 'z';
    this->k_toggle_experimental_ = 'u';
    this->k_inc_thickness_ = '+';
    this->k_dec_thickness_ = '-';
    this->k_save_skel_ = -2;
    this->k_save_last_poly_ = -3;
    this->k_dump_window_ = -10;
    this->k_print_screen_ = -11;
    this->k_print_cutpattern_ = -12;

    this->speed_move_ = 0.1;
    this->speed_look_ = 0.02;

    loadConfig(util::Configuration::getInstance());

    this->move_forward_ = false;
    this->move_backpedal_ = false;
    this->move_left_ = false;
    this->move_right_ = false;
    this->move_up_ = false;
    this->move_down_ = false;
    this->look_left_ = false;
    this->look_right_ = false;
    this->look_up_ = false;
    this->look_down_ = false;

    this->thread_ = ThreadSPtr(new std::thread(
            std::bind(&KeyboardAdapter::run, this)));
}

KeyboardAdapter::~KeyboardAdapter() {
    this->camera_.reset();
}

KeyboardAdapterSPtr KeyboardAdapter::create(CameraSPtr camera, ControllerSPtr controller) {
    KeyboardAdapterSPtr result = KeyboardAdapterSPtr(new KeyboardAdapter(camera, controller));
    return result;
}

void KeyboardAdapter::setWindow(MainOpenGLWindowSPtr window) {
    this->window_ = window;
}

int KeyboardAdapter::toKey(const std::string& key) {
    int result = 0;
    if (key.length() == 1) {
        result = (int)key.c_str()[0];
    } else if (key.compare("ESC") == 0) {
        result = (int)'\033';
    } else if (key.compare("SPACE") == 0) {
        result = (int)' ';
    } else if (key.compare("BACKSPACE") == 0) {
        result = 8;
    } else if (key.compare("UPARROW") == 0) {
        result = -101;
    } else if (key.compare("DOWNARROW") == 0) {
        result = -103;
    } else if (key.compare("LEFTARROW") == 0) {
        result = -100;
    } else if (key.compare("RIGHTARROW") == 0) {
        result = -102;
    } else if (key.compare("F1") == 0) {
        result = -1;
    } else if (key.compare("F2") == 0) {
        result = -2;
    } else if (key.compare("F3") == 0) {
        result = -3;
    } else if (key.compare("F4") == 0) {
        result = -4;
    } else if (key.compare("F5") == 0) {
        result = -5;
    } else if (key.compare("F6") == 0) {
        result = -6;
    } else if (key.compare("F7") == 0) {
        result = -7;
    } else if (key.compare("F8") == 0) {
        result = -8;
    } else if (key.compare("F9") == 0) {
        result = -9;
    } else if (key.compare("F10") == 0) {
        result = -10;
    } else if (key.compare("F11") == 0) {
        result = -11;
    } else if (key.compare("F12") == 0) {
        result = -12;
    }
    return result;
}

void KeyboardAdapter::loadConfig(util::ConfigurationSPtr config) {
    if (!config->isLoaded()) {
        return;
    }
    std::string section("ui_gl_KeyboardAdapter");
    int key = 0;
    key = toKey(config->getString(section, "k_move_forward"));
    if (key != 0) k_move_forward_ = key;
    key = toKey(config->getString(section, "k_move_backpedal"));
    if (key != 0) k_move_backpedal_ = key;
    key = toKey(config->getString(section, "k_move_left"));
    if (key != 0) k_move_left_ = key;
    key = toKey(config->getString(section, "k_move_right"));
    if (key != 0) k_move_right_ = key;
    key = toKey(config->getString(section, "k_move_up"));
    if (key != 0) k_move_up_ = key;
    key = toKey(config->getString(section, "k_move_down"));
    if (key != 0) k_move_down_ = key;
    key = toKey(config->getString(section, "k_look_left"));
    if (key != 0) k_look_left_ = key;
    key = toKey(config->getString(section, "k_look_right"));
    if (key != 0) k_look_right_ = key;
    key = toKey(config->getString(section, "k_look_up"));
    if (key != 0) k_look_up_ = key;
    key = toKey(config->getString(section, "k_look_down"));
    if (key != 0) k_look_down_ = key;
    key = toKey(config->getString(section, "k_reset"));
    if (key != 0) k_reset_ = key;
    key = toKey(config->getString(section, "k_pause"));
    if (key != 0) k_pause_ = key;
    key = toKey(config->getString(section, "k_step"));
    if (key != 0) k_step_ = key;
    key = toKey(config->getString(section, "k_skip"));
    if (key != 0) k_skip_ = key;
    key = toKey(config->getString(section, "k_toggle_roof"));
    if (key != 0) k_toggle_roof_ = key;
    key = toKey(config->getString(section, "k_toggle_poly"));
    if (key != 0) k_toggle_poly_ = key;
    key = toKey(config->getString(section, "k_toggle_skel"));
    if (key != 0) k_toggle_skel_ = key;
    key = toKey(config->getString(section, "k_toggle_experimental"));
    if (key != 0) k_toggle_experimental_ = key;
    key = toKey(config->getString(section, "k_inc_thickness"));
    if (key != 0) k_inc_thickness_ = key;
    key = toKey(config->getString(section, "k_dec_thickness"));
    if (key != 0) k_dec_thickness_ = key;
    key = toKey(config->getString(section, "k_save_skel"));
    if (key != 0) k_save_skel_ = key;
    key = toKey(config->getString(section, "k_save_last_poly"));
    if (key != 0) k_save_last_poly_ = key;
    key = toKey(config->getString(section, "k_dump_window"));
    if (key != 0) k_dump_window_ = key;
    key = toKey(config->getString(section, "k_print_screen"));
    if (key != 0) k_print_screen_ = key;
    key = toKey(config->getString(section, "k_print_cutpattern"));
    if (key != 0) k_print_cutpattern_ = key;

    double speed = 0.0;
    speed = config->getDouble(section, "speed_move");
    if (speed != 0.0) speed_move_ = speed;
    speed = config->getDouble(section, "speed_look");
    if (speed != 0.0) speed_look_ = speed;
}

void KeyboardAdapter::pressed(int key, int x, int y) {
    // DEBUG_VAR(key);
    if (key == k_move_forward_) move_forward_ = true;
    if (key == k_move_backpedal_) move_backpedal_ = true;
    if (key == k_move_left_) move_left_ = true;
    if (key == k_move_right_) move_right_ = true;
    if (key == k_move_down_) move_down_ = true;
    if (key == k_move_up_) move_up_ = true;
    if (key == k_look_left_) look_left_ = true;
    if (key == k_look_right_) look_right_ = true;
    if (key == k_look_up_) look_up_ = true;
    if (key == k_look_down_) look_down_ = true;
    if (key == k_reset_) this->camera_->reset();
    if (key == k_pause_) this->controller_->togglePause();
    if (key == k_step_) this->controller_->nextStep();
    if (key == k_skip_) this->controller_->skip();
    if (key == k_toggle_roof_) {
        if (!this->window_.expired()) {
            this->window_.lock()->toggleRoof();
            this->camera_->topdown();
        }
    }
    if (key == k_toggle_poly_) {
        if (!this->window_.expired()) this->window_.lock()->togglePoly();
    }
    if (key == k_toggle_skel_) {
        if (!this->window_.expired()) this->window_.lock()->toggleSkel();
    }
    if (key == k_toggle_experimental_) {
        if (!this->window_.expired()) this->window_.lock()->toggleExperimental();
    }
    if (key == k_inc_thickness_) {
        if (!this->window_.expired()) this->window_.lock()->incThickness();
    }
    if (key == k_dec_thickness_) {
        if (!this->window_.expired()) this->window_.lock()->decThickness();
    }
    if (key == k_save_skel_) {
        if (!this->window_.expired()) this->window_.lock()->saveSkel();
    }
    if (key == k_save_last_poly_) {
        if (!this->window_.expired()) this->window_.lock()->saveLastPoly();
    }
    if (key == k_dump_window_) {
        if (!this->window_.expired()) this->window_.lock()->dumpWin();
    }
    if (key == k_print_screen_) {
        if (!this->window_.expired()) this->window_.lock()->printScreen();
    }
    if (key == k_print_cutpattern_) {
        if (!this->window_.expired()) this->window_.lock()->printCutPattern();
    }
}

void KeyboardAdapter::released(int key, int x, int y) {
    if (key == k_move_forward_) move_forward_ = false;
    if (key == k_move_backpedal_) move_backpedal_ = false;
    if (key == k_move_left_) move_left_ = false;
    if (key == k_move_right_) move_right_ = false;
    if (key == k_move_down_) move_down_ = false;
    if (key == k_move_up_) move_up_ = false;
    if (key == k_look_left_) look_left_ = false;
    if (key == k_look_right_) look_right_ = false;
    if (key == k_look_up_) look_up_ = false;
    if (key == k_look_down_) look_down_ = false;
}

void KeyboardAdapter::run() {
    while (true) {
        if (move_forward_) camera_->moveFB(speed_move_);
        if (move_backpedal_) camera_->moveFB(-speed_move_);
        if (move_left_) camera_->strafeLR(-speed_move_);
        if (move_right_) camera_->strafeLR(speed_move_);
        if (move_up_) camera_->moveUD(speed_move_);
        if (move_down_) camera_->moveUD(-speed_move_);
        if (look_left_) camera_->lookLR(-speed_look_);
        if (look_right_) camera_->lookLR(speed_look_);
        if (look_up_) camera_->lookUD(speed_look_);
        if (look_down_) camera_->lookUD(-speed_look_);
        thread_sleep(10);
    }
}

} }
