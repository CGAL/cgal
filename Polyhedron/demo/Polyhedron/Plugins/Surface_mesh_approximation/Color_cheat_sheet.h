#ifndef COLOR_CHEAT_SHEET_H
#define COLOR_CHEAT_SHEET_H

#include <cstddef>

// 256 color cheat sheet
// Source: https://jonasjacek.github.io/colors/

#ifdef color_cheat_sheet_EXPORTS
# define COLOR_CHEAT_SHEET_EXPORT Q_DECL_EXPORT
#else
# define COLOR_CHEAT_SHEET_EXPORT Q_DECL_IMPORT
#endif

class COLOR_CHEAT_SHEET_EXPORT Color_cheat_sheet
{
public:
  static const unsigned char &r(const std::size_t &i) {
    return m_colors[i][0];
  }

  static const unsigned char &g(const std::size_t &i) {
    return m_colors[i][1];
  }

  static const unsigned char &b(const std::size_t &i) {
    return m_colors[i][2];
  }

private:
  static const unsigned char m_colors[256][3];
};

#endif // COLOR_CHEAT_SHEET_H
