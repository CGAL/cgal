#version 150

in GS_OUT
{
 vec4 color;
} fs_in;

uniform float width;

out vec4 out_color;

void main(void)
{
  out_color = fs_in.color;
}
