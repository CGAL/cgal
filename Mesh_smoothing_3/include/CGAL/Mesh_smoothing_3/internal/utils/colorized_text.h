#pragma once

// taken from https://stackoverflow.com/a/78739378

#include <cstdint>
#include <iostream>
#if defined WIN32 || defined _WIN64
#include <Windows.h>
#endif

namespace Mesh_smoothing_3_internal {

#if defined WIN32 || defined _WIN64
enum class ConsoleTextColor : std::uint8_t
{
    Default = 7,
    Blue = 1,
    Green = 2,
    Cyan = 3,
    Red = 4,
    Pink = 5,
    Yellow = 6,
    White = 7,
    BrightBlue = 9,
    BrightGreen = 10,
    BrightCyan = 11,
    BrightRed = 12,
    BrightPink = 13,
    BrightYellow = 14,
    BrightWhite = 15
};
#else

enum class ConsoleTextColor : std::uint8_t
{
    Default = 33,
    Blue = 34,
    Green = 32,
    Cyan = 36,
    Red = 31,
    Pink = 35,
    Yellow = 33,
    White = 37,
    BrightBlue = 94,
    BrightGreen = 92,
    BrightCyan = 96,
    BrightRed = 91,
    BrightPink = 95,
    BrightYellow = 93,
    BrightWhite = 97
};

#endif

inline void Colorized_print(const std::string &text, const ConsoleTextColor& color = ConsoleTextColor::Default) {

    if (color != ConsoleTextColor::Default)
    {
        #if defined WIN32 || defined _WIN64

        // Get console descriptor
        HANDLE ConsoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);

        // Object to store previous console colour
        CONSOLE_SCREEN_BUFFER_INFO info;
        // Tty to store console colour
        if (GetConsoleScreenBufferInfo(ConsoleHandle, &info))
        {
            SetConsoleTextAttribute(ConsoleHandle, static_cast<std::uint8_t>(color));
            std::cout << text << std::endl;
            SetConsoleTextAttribute(ConsoleHandle, static_cast<std::uint8_t>(info.wAttributes));
        }
        #else
        std::cout << "\x1B[" + std::to_string(static_cast<std::uint8_t>(color)) + "m" + text + "\033[0m" << std::endl;
        #endif
    }
    else
        std::cout << text << std::endl;

}

}
