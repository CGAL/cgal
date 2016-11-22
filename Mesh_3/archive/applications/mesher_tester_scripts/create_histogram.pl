#!/usr/bin/perl

BEGIN {
  push(@INC,$ENV{'HOME'}."/local/perl"); #where the SVG.pm file lives
  push(@INC,$ENV{'HOME'}."/local/perl/SVG"); # where the Utils.pm file lives
}


use strict;
use SVG;

my $height = 140;
my $width = 220;
my $shift = 20;
my $padding_v = 5;
my $svg=new SVG(width=>$width,height=>$height+30, style=>{'background'=>'lightgrey'});

my $file = @ARGV[0];
if ( -f $file ) {
  open(my $in, $file) || die("Impossible d'ouvrir $file: $!");

  my @values;
  my $max_value = 0;
  while (<$in>) {
    my $line = $_;
    chomp($line);
    if ( $line > $max_value ) {
      $max_value = $line;
    }

    push(@values,$line);
  }


  # normalize values
  my $min = shift(@values);
  my $max = shift(@values);
  my @norm = map( { $_/$max_value*$height } @values);


  # draw
  my $i = 0;
  my $color = "black";
  foreach my $v (@norm) {
    # use red for the isotrop tet
    if ( $i == 71 ) {
      $color = "red";
    }
    
    my $tag = $svg->rectangle(x=>$i+$shift, y=>$height-$v+$padding_v,
                              width=>1, height=>$v,
                              rx=>0, ry=>0,
                              id=>"rect_$i",
                              style=>{'fill'=>$color});
    $color = "black";
    $i++;
  }
  
  $svg->text(id=>'min', x=>5, y=>$height+20)->cdata($min);
  $svg->text(id=>'max', x=>160, y=>$height+20)->cdata($max);
}
else {
  $svg->text(id=>'error', x=>$width/2-50, y=>$height/2+20, style=>{'fill'=>'black'})->cdata("No data found...");
}

# now render the SVG object, implicitly use svg namespace
print $svg->xmlify;
