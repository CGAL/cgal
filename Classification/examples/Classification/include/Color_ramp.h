#ifndef COLOR_RAMP_H
#define COLOR_RAMP_H

class Color_ramp
{
  typedef std::array<unsigned char, 3> Color;
  typedef std::pair<double, Color> Step;

  std::vector<Step> m_steps;

public:

  Color_ramp()
  {
    m_steps.push_back (std::make_pair (0, Color{ 192, 192, 255}));
    m_steps.push_back (std::make_pair (0.2, Color{ 0, 0, 255}));
    m_steps.push_back (std::make_pair (0.4, Color{ 0, 255, 0}));
    m_steps.push_back (std::make_pair (0.6, Color{ 255, 255, 0}));
    m_steps.push_back (std::make_pair (0.8, Color{ 255, 0, 0}));
    m_steps.push_back (std::make_pair (1.0, Color{ 128, 0, 0}));
  }

  Color get (double value) const
  {
    std::size_t idx = 0;
    while (m_steps[idx+1].first < value)
      ++ idx;

    double v0 = m_steps[idx].first;
    double v1 = m_steps[idx+1].first;
    const Color& c0 = m_steps[idx].second;
    const Color& c1 = m_steps[idx+1].second;

    double ratio = (value - v0) / (v1 - v0);

    Color out;
    for (std::size_t i = 0; i < 3; ++ i)
      out[i] = (unsigned char)((1 - ratio) * c0[i] + ratio * c1[i]);

    return out;
  }

};

#endif // COLOR_RAMP_H
