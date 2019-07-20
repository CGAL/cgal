#version 150
in vec3 fColors;
in vec2 f_texCoord; 
uniform sampler2D s_texture; 
out vec4 out_color;
 
void main(void) 
{ 
  out_color = vec4(vec3(texture(s_texture, f_texCoord))*fColors, 1.0);
}  
