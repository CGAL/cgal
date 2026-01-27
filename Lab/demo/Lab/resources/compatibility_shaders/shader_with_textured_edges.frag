
varying highp vec3 fColors;
varying highp vec2 f_texCoord; 
uniform highp sampler2D s_texture; 
 
void main(void) 
{ 
  gl_FragColor = vec4(vec3(texture2D(s_texture, f_texCoord))*fColors, 1.0); 
}
