#version 330
varying highp vec3 fColors;
varying highp vec2 f_texCoord; 
uniform sampler2D s_texture; 
 
void main(void) 
{ 
  gl_FragColor = vec4(vec3(texture(s_texture, f_texCoord))*fColors, 1.0);
}  
