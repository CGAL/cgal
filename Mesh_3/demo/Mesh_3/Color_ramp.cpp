#include "Color_ramp.h"


Color_ramp::Color_ramp()
{ 
  build_thermal(); 
}

void 
Color_ramp::rebuild()
{
  // build nodes
  m_colors[3][0] = 1;
  m_colors[3][255] = 1;
  unsigned int nb_nodes = 0;
  for(int i=0;i<256;i++)
  {
    if(m_colors[3][i])
    {
      m_nodes[nb_nodes] = i;
      nb_nodes++;
    }
  }
  
  // build color_ramp
  for(int k=0;k<3;k++)
    for(unsigned int i=0;i<(nb_nodes-1);i++)
    {
      int x1 = m_nodes[i];
      int x2 = m_nodes[i+1];
      int y1 = m_colors[k][x1];
      int y2 = m_colors[k][x2];
      float a = (float)(y2-y1) / (float)(x2-x1);
      float b = (float)y1 - a*(float)x1;
      for(int j=x1;j<x2;j++)
        m_colors[k][j] = (unsigned char)(a*(float)j+b);
    }
}

void 
Color_ramp::reset()
{
  for(int i=1;i<=254;i++)
    m_colors[3][i] = 0;
  m_colors[3][0] = 1;
  m_colors[3][255] = 1;
}

void 
Color_ramp::add_node(unsigned int index,
                     unsigned char r,
                     unsigned char g,
                     unsigned char b)
{
  m_colors[3][index] = 1;
  m_colors[0][index] = r;
  m_colors[1][index] = g;   
  m_colors[2][index] = b;
}

void 
Color_ramp::build_red()
{
  reset();
  m_colors[3][0] = 1;
  m_colors[0][0] = 128;
  m_colors[1][0] = 0;
  m_colors[2][0] = 0;
  
  m_colors[3][80] = 1;
  m_colors[0][80] = 255;
  m_colors[1][80] = 0;
  m_colors[2][80] = 0;
  
  m_colors[3][160] = 1;
  m_colors[0][160] = 255;
  m_colors[1][160] = 128;
  m_colors[2][160] = 0;
  
  m_colors[3][255] = 1;
  m_colors[0][255] = 255;
  m_colors[1][255] = 255;
  m_colors[2][255] = 255;
  
  rebuild();
}

void 
Color_ramp::build_thermal()
{
  reset();
  m_colors[3][0] = 1;
  m_colors[0][0] = 0;
  m_colors[1][0] = 0;
  m_colors[2][0] = 0;
  
  m_colors[3][70] = 1;
  m_colors[0][70] = 128;
  m_colors[1][70] = 0;
  m_colors[2][70] = 0;
  
  m_colors[3][165] = 1;
  m_colors[0][165] = 255;
  m_colors[1][165] = 128;
  m_colors[2][165] = 0;
  
  m_colors[3][240] = 1;
  m_colors[0][240] = 255;
  m_colors[1][240] = 191;
  m_colors[2][240] = 128;
  
  m_colors[3][255] = 1;
  m_colors[0][255] = 255;
  m_colors[1][255] = 255;
  m_colors[2][255] = 255;
  
  rebuild();
}

void 
Color_ramp::build_blue()
{
  reset();
  m_colors[3][0] = 1;
  m_colors[0][0] = 0;
  m_colors[1][0] = 0;
  m_colors[2][0] = 128;
  
  m_colors[3][80] = 1;
  m_colors[0][80] = 0;
  m_colors[1][80] = 0;
  m_colors[2][80] = 255;
  
  m_colors[3][160] = 1;
  m_colors[0][160] = 0;
  m_colors[1][160] = 255;
  m_colors[2][160] = 255;
  
  m_colors[3][255] = 1;
  m_colors[0][255] = 255;
  m_colors[1][255] = 255;
  m_colors[2][255] = 255;
  
  rebuild();
}
