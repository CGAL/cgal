#version 150

in vec4 vertex;
in vec4 colors;

uniform mat4 mvp_matrix;

out VS_OUT
{
  vec4 out_color;
}vs_out;

void main(void)
{
   gl_Position = mvp_matrix * vertex;
   vs_out.out_color = colors;
}
