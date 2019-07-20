#version 150

in GS_OUT
{
 vec4 color;
 float dist[6];
} fs_in;

uniform float width;
uniform bool is_selected;
uniform bool is_clipbox_on;


out vec4 out_color;

void main(void)
{
  if(is_clipbox_on)
  {
    if(fs_in.dist[0]>0.0 ||
       fs_in.dist[1]>0.0 ||
       fs_in.dist[2]>0.0 ||
       fs_in.dist[3]>0.0 ||
       fs_in.dist[4]>0.0 ||
       fs_in.dist[5]>0.0)
      discard;
    else
      out_color = fs_in.color;
   }
  else
    out_color = fs_in.color;
}
