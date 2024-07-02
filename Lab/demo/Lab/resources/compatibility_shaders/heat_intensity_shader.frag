
varying highp vec4 color;
varying highp float out_dist;
void main(void) {
  highp vec3 c = color.xyz;
  highp float h = out_dist;
  h = h * 20.;
  h = h - floor(h);
  h = (1./(1.+exp(-100.*(h-.55)))) + (1./(1.+exp(-100.*(-h+.45))));
  h = 1.-h;

  c = h*vec3(1.,1.,1.) + (1.-h)*c;
  
  gl_FragColor.rgb = c;
  gl_FragColor.a = 1.;
}
