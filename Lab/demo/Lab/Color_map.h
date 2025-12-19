#ifndef _COLOR_MAP_H
#define _COLOR_MAP_H

#include <QColor>

#include <stdlib.h>
#include <iostream>

inline QColor generate_color(double h,
                             double s_min = 0.35)
{
  std::size_t s_max = 255;
  if(h > 0.8 && h < 0.95) // span of ugly pink, desaturates make it less garish IMO
    s_max = 160;
  std::size_t s = std::rand() % (s_max-static_cast<std::size_t>(s_min*255)) + static_cast<int>(s_min*255);
  return QColor::fromHsvF(h, s/255.0, 1.0);
}


template <typename Output_color_iterator>
Output_color_iterator
compute_color_map(QColor base_color,
                  std::size_t nb_of_colors,
                  Output_color_iterator out)
{
  const qreal step = (static_cast<qreal>(0.85)) / nb_of_colors;

  qreal hue = base_color.hueF();
  qreal h = (hue == -1) ? 0 : hue;
  for(std::size_t i=0; i<nb_of_colors; ++i)
  {
    if(h != -1)
      h += step;
    if(h > 1)
      h -= 1;
    *out++ = generate_color(h);
  }

  return out;
}

template <typename Output_color_iterator>
Output_color_iterator
compute_deterministic_color_map(QColor base_color,
                                std::size_t nb_of_colors,
                                Output_color_iterator out)
{
  qreal hue = base_color.hueF();
  qreal saturation = base_color.saturationF();
  qreal value = base_color.valueF();
  const qreal hue_step = (hue == -1) ? 0 : (static_cast<qreal>(1)) / nb_of_colors;

  if (hue == -1)
    hue = 0;
  for(std::size_t i=0; i<nb_of_colors; ++i)
  {
    hue += hue_step;
    if(hue > 1)
      hue -= 1;
    *out++ = QColor::fromHsvF(hue, saturation, value);
  }

  return out;
}

inline QColor generate_random_color()
{
  std::size_t h = static_cast<std::size_t>(std::rand() % 360);
  return generate_color(h / 359.0);
}

#endif // _COLOR_MAP_H
