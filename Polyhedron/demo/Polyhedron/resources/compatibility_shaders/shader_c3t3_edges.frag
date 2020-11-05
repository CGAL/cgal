
varying highp vec4 color;
void main(void)
{
  if(color.w<0.)
    gl_FragColor = vec4(0.,0.,0.,1.0);
  else
    discard;
}
