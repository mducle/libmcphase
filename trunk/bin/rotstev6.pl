#!/usr/bin/perl
# rotstev.pl 
#
# Generates rotation matrices to rotate Stevens' Operator equivalents by
# phi about y-axis and theta about z-axis. Based on the method of Buckmaster
# (Physica Status Solidi A vol 13, pp 9, 1972) and Rudowicz (J. Phys: C: 
# Solid State Physics, vol 18, pp 1415, 1985).
#
# NB: phi is denoted by f and theta by t in this program
#
# by Duc Le 2006 - duc.le@ucl.ac.uk

use Math::Algebra::Symbols trig=>1;     # exports trig function to my namespace
use Getopt::Long;

# The following subroutines were shamelessly stolen from the Perl Cookbook
# by O'Reily, and used for the symbolic matrix multiplications.

sub mmult {
  my ($m1,$m2) = @_;
  my ($m1rows,$m1cols) = matdim($m1);
  my ($m2rows,$m2cols) = matdim($m2);

  unless ($m1cols == $m2rows) {       # raise exception
    die "IndexError: matrices don't match: $m1cols != $m2rows";
  }

  my $result = [];
  my ($i, $j, $k, $m1el, $m2el, $prod_el);

  for $i (range($m1rows)) {
    for $j (range($m2cols)) {
      for $k (range($m1cols)) {
        $m1el = $m1->[$i][$k];
        $m2el = $m2->[$k][$j];
        $prod_el += $m1el * $m2el;
      }
      $result->[$i][$j] = $prod_el;
      $prod_el = 0;
    }
  }

  return $result;
}

sub range { 0 .. ($_[0] - 1) }

sub veclen {
    my $ary_ref = $_[0];
    my $type = ref $ary_ref;
    if ($type ne "ARRAY") { die "$type is bad array ref for $ary_ref" }
    return scalar(@$ary_ref);
}

sub matdim {
    my $matrix = $_[0];
    my $rows = veclen($matrix);
    my $cols = veclen($matrix->[0]);
    return ($rows, $cols);
}

# End of stolen code!

sub usage() {
  print STDERR << "EOF";

    $0:
    Calculates the relationship  between a set of crystal field  parameters
    which has  been rotated by phi about  the y-axis and theta about the z-
    axis and its unrotated counter-part.

    The  calculations  are  done  by  means  of  matrix  multiplication and 
    symbolic algebra  using the Math::Algebra::Symbols package and based on
    the method of  Buckmaster (phys. stat. sol. a, vol 13,  pp 9, 1972) and
    Rudowicz (J. Phys: Solid State Phys., vol 18, pp 1415, 1985).   

    usage: $0 [-h] [--help] 
              [-r] [-radians]
	      [-p phi_angle] [--phi phi_angle]
	      [-t theta_angle] [--theta theta_angle]
	      [-o output_file] [--output output_file]

     -h             : this (help) message
     -r             : specifies that phi and theta are given in radians
     -p phi_angle   : angle to rotate about y-axis (=0 by default)
     -t theta_angel : angle to rotate about z-axis (=90 degrees by default)  
     -o out_file    : output latex file name

    if -r is omitted, the program assumes phi and theta are in degrees.
    
    if either -p or -t are omitted, the program will use the default values 
    
    if -o is omitted the  program will output the  transformation relations
    and rotation matrix to screen.  Otherwise it will generate a Latex file
    of the transformations and matrix.

    (C) 2006 duc.le\@ucl.ac.uk
EOF
  exit;
}

# see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions

GetOptions("help"=>\$helpflag,
	   "radians"=>\$radflag,
	   "phi=f"=>\$phi,
	   "theta=f"=>\$theta,
	   "output=s"=>\$output,
           "z_internal"=>\$z_int);
	   
usage() if $helpflag;

my ($S, $f, $t, $i, $pi) = symbols(qw(six f t i pi));
my ($l, $m, $n);

if (!$phi) { $phi = 0; }
if (!$theta) { $theta = 90; }

if (!$radflag) {
  $phi_val = $phi * $pi / 180;
  $theta_val = $theta * $pi / 180;
}
else {
  $phi_val = $phi;
  $theta_val = $theta;
}

my ($f, $t, $i, $pi) = symbols(qw(f t i pi));
my ($T, $e, $F, $X, $S, $E) = symbols(qw(two three five six seven eleven));
my ($l, $m, $n);

my $A6 = [
	[4,    0,     0,    0,     0,     0,     0,      0,      0,     0,     0,     0,     -4],
	[0,2/sqrt($e),0,    0,     0,     0,     0,      0,      0,     0,     0, 2/sqrt($e), 0],
	[0,    0,4*sqrt($E/$X),0,  0,     0,     0,      0,      0,     0,-4*sqrt($E/$X),0,   0],
	[0,    0,     0,2*sqrt($E/$F),0,  0,     0,      0,      0,2*sqrt($E/$F),0,   0,      0],
	[0,    0,     0,    0,4*sqrt($E/$F),0,   0,      0,-4*sqrt($E/$F),0,   0,     0,      0],
	[0,    0,     0,    0,     0,sqrt($T*$E),0,  sqrt($T*$E),0,     0,     0,     0,      0],
        [0,    0,     0,    0,     0,   0,4*sqrt($e*$S*$E),0,    0,     0,     0,     0,      0],
	[0,    0,     0,    0,     0,sqrt($T*$E),0, -sqrt($T*$E),0,     0,     0,     0,      0],
	[0,    0,     0,    0,4*sqrt($E/$F),0,   0,      0, 4*sqrt($E/$F),0,   0,     0,      0],
	[0,    0,     0,2*sqrt($E/$F),0,  0,     0,      0,      0,-2*sqrt($E/$F),0,  0,      0],
	[0,    0,4*sqrt($E/$X),0,  0,     0,     0,      0,      0,     0, 4*sqrt($E/$X),0,   0],
	[0,2/sqrt($e),0,    0,     0,     0,     0,      0,      0,     0,     0,-2/sqrt($e), 0],
	[4,    0,     0,    0,     0,     0,     0,      0,      0,     0,     0,     0,      4]
      ];

my $invA6 = [
	[  1/8,                0,0,0,0,0,0,0,0,0,0,0,      1/8           ],
	[0,  sqrt($e)/4,         0,0,0,0,0,0,0,0,0,    sqrt($e)/4,      0],
	[0,0,  sqrt($e/$T/$E)/4,   0,0,0,0,0,0,0,   sqrt($e/$T/$E)/4, 0,0],
	[0,0,0,  sqrt($F/$E)/4,      0,0,0,0,0,   sqrt($F/$E)/4,    0,0,0],
	[0,0,0,0,  sqrt($F/$E)/8,     0,0,0,    sqrt($F/$E)/8,    0,0,0,0],
	[0,0,0,0,0,  sqrt(1/$T/$E)/2,   0,    sqrt(1/$T/$E)/2,  0,0,0,0,0],
	[0,0,0,0,0,0,         sqrt(1/$e/$S/$E)/4,             0,0,0,0,0,0],
	[0,0,0,0,0,   sqrt(1/$T/$E)/2,  0,   -sqrt(1/$T/$E)/2,  0,0,0,0,0],
	[0,0,0,0, -sqrt($F/$E)/8,     0,0,0,   sqrt($F/$E)/8,     0,0,0,0],
	[0,0,0,  sqrt($F/$E)/4,     0,0,0,0,0,  -sqrt($F/$E)/4,     0,0,0],
	[0,0, -sqrt($e/$T/$E)/4,  0,0,0,0,0,0,0,   sqrt($e/$T/$E)/4,  0,0],
	[0,  sqrt($e)/4,        0,0,0,0,0,0,0,0,0,  -sqrt($e)/4,        0],
	[ -1/8,               0,0,0,0,0,0,0,0,0,0,0,       1/8            ]
      ];

# --------------------------------------- O6 ------------------------------------- #

my $D60p0 =  (1/16)            * (231*cos($t)**6 - 315*cos($t)**4 + 105*cos($t)**2 - 5);
my $D60p1 =  sqrt($X*$S)/16    * sin($t)   * cos($t) * (33*cos($t)**4 - 30*cos($t)**2 + 5);
my $D60p2 =  sqrt($e*$F*$S)/32 * sin($t)**2 * (33*cos($t)**4 - 18*cos($t)**2 + 1);
my $D60p3 =  sqrt($e*$F*$S)/16 * sin($t)**3 * cos($t) * (11*cos($t)**2 - 3);
my $D60p4 =  3*sqrt($T*$S)/32  * sin($t)**4 * (11*cos($t)**2 - 1);
my $D60p5 =  3*sqrt($S*$E)/16  * sin($t)**5 * cos($t);
my $D60p6 =  sqrt($e*$S*$E)/32 * sin($t)**6;

my $D61p0 =  exp( $i*$f) * sqrt($X*$S)/16   * sin($t)   * cos($t) * (33*cos($t)**4 - 30*cos($t)**2 + 5);
my $D61p1 =  exp( $i*$f) * (1/16) * (1+cos($t)) * (198*cos($t)**5 - 165*cos($t)**4 - 120*cos($t)**3 + 90*cos($t)**2 + 10*cos($t) - 5);
my $D61m1 =  exp(-$i*$f) * (1/16) * (1-cos($t)) * (198*cos($t)**5 + 165*cos($t)**4 - 120*cos($t)**3 - 90*cos($t)**2 + 10*cos($t) + 5);
my $D61p2 =  exp( $i*$f) * sqrt($T*$F)/32   * sin($t)   * (1+cos($t)) * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1);
my $D61m2 = -exp(-$i*$f) * sqrt($T*$F)/32   * sin($t)   * (1-cos($t)) * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1);
my $D61p3 =  exp( $i*$f) * 3*sqrt($T*$F)/32 * sin($t)**2 * (1+cos($t)) * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D61m3 =  exp(-$i*$f) * 3*sqrt($T*$F)/32 * sin($t)**2 * (1-cos($t)) * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D61p4 =  exp( $i*$f) * sqrt($e)/16      * sin($t)**3 * (1+cos($t)) * (33*cos($t)**2 - 11*cos($t) - 2);
my $D61m4 = -exp(-$i*$f) * sqrt($e)/16      * sin($t)**3 * (1-cos($t)) * (33*cos($t)**2 + 11*cos($t) - 2);
my $D61p5 =  exp( $i*$f) * sqrt($X*$E)/32   * sin($t)**4 * (1+cos($t)) * (6*cos($t) - 1);
my $D61m5 =  exp(-$i*$f) * sqrt($X*$E)/32   * sin($t)**4 * (1-cos($t)) * (6*cos($t) + 1);
my $D61p6 =  exp( $i*$f) * 3*sqrt($T*$E)/32 * sin($t)**5 * (1+cos($t));
my $D61m6 = -exp(-$i*$f) * 3*sqrt($T*$E)/32 * sin($t)**5 * (1-cos($t));

my $D62p0 =  exp( 2*$i*$f) * sqrt($e*$F*$S)/32 * sin($t)**2 * (33*cos($t)**4 - 18*cos($t)**2 + 1);
my $D62p1 =  exp( 2*$i*$f) * sqrt($T*$F)/32    * sin($t)    * (1+cos($t))   * (99*cos($t)**4 - 66*cos($t)**3 - 36*cos($t)**2 + 18*cos($t) + 1); 
my $D62m1 =  exp(-2*$i*$f) * sqrt($T*$F)/32    * sin($t)    * (1-cos($t))   * (99*cos($t)**4 + 66*cos($t)**3 - 36*cos($t)**2 - 18*cos($t) + 1); 
my $D62p2 =  exp( 2*$i*$f) * (1/64) *  (1+cos($t))**2 * (495*cos($t)**4 - 660*cos($t)**3 + 90*cos($t)**2 + 108*cos($t) - 17);
my $D62m2 =  exp(-2*$i*$f) * (1/64) *  (1-cos($t))**2 * (495*cos($t)**4 + 660*cos($t)**3 + 90*cos($t)**2 - 108*cos($t) - 17);
my $D62p3 =  exp( 2*$i*$f) * (3/32)            * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D62m3 = -exp(-2*$i*$f) * (3/32)            * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D62p4 =  exp( 2*$i*$f) * sqrt($X*$F)/64    * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D62m4 =  exp(-2*$i*$f) * sqrt($X*$F)/64    * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D62p5 =  exp( 2*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D62m5 = -exp(-2*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D62p6 =  exp( 2*$i*$f) * 3*sqrt($F*$E)/64  * sin($t)**4 * (1+cos($t))**2;
my $D62m6 =  exp(-2*$i*$f) * 3*sqrt($F*$E)/64  * sin($t)**4 * (1-cos($t))**2;

my $D63p0 =  exp( 3*$i*$f) * sqrt($e*$F*$E)/16 * sin($t)**3 * cos($t)        * (11*cos($t)**2 - 3);
my $D63p1 =  exp( 3*$i*$f) * 3*sqrt($T*$F)/32  * sin($t)**2 * (1+cos($t))    * (22*cos($t)**3 - 11*cos($t)**2 - 4*cos($t) + 1);
my $D63m1 =  exp(-3*$i*$f) * 3*sqrt($T*$F)/32  * sin($t)**2 * (1-cos($t))    * (22*cos($t)**3 + 11*cos($t)**2 - 4*cos($t) - 1);
my $D63p2 =  exp( 3*$i*$f) * (3/32)            * sin($t)    * (1+cos($t))**2 * (55*cos($t)**3 - 55*cos($t)**2 + 5*cos($t) + 3);
my $D63m2 =  exp(-3*$i*$f) * (3/32)            * sin($t)    * (1-cos($t))**2 * (55*cos($t)**3 + 55*cos($t)**2 + 5*cos($t) - 3);
my $D63p3 =  exp( 3*$i*$f) * (1/32)            *           (1+cos($t))**3 * (110*cos($t)**3 - 165*cos($t)**2 + 60*cos($t) - 1);
my $D63m3 =  exp(-3*$i*$f) * (1/32)            *           (1-cos($t))**3 * (110*cos($t)**3 + 165*cos($t)**2 + 60*cos($t) + 1);
my $D63p4 =  exp( 3*$i*$f) * sqrt($F*$X)/32    * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D63m4 = -exp(-3*$i*$f) * sqrt($F*$X)/32    * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D63p5 =  exp( 3*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D63m5 =  exp(-3*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D63p6 =  exp( 3*$i*$f) * sqrt($F*$E)/32    * sin($t)**3 * (1+cos($t))**3;
my $D63m6 = -exp(-3*$i*$f) * sqrt($F*$E)/32    * sin($t)**3 * (1-cos($t))**3;
 
my $D64p0 =  exp( 4*$i*$f) * 3*sqrt($T*$S)/32 * sin($t)**4 * (11*cos($t)**2 - 1);
my $D64p1 =  exp( 4*$i*$f) * sqrt($e)/16      * sin($t)**3 * (1+cos($t))    * (33*cos($t)**2 - 11*cos($t) - 2);
my $D64m1 =  exp(-4*$i*$f) * sqrt($e)/16      * sin($t)**3 * (1-cos($t))    * (33*cos($t)**2 + 11*cos($t) - 2);
my $D64p2 =  exp( 4*$i*$f) * sqrt($X*$F)/64   * sin($t)**2 * (1+cos($t))**2 * (33*cos($t)**2 - 22*cos($t) + 1);
my $D64m2 =  exp(-4*$i*$f) * sqrt($X*$F)/64   * sin($t)**2 * (1-cos($t))**2 * (33*cos($t)**2 + 22*cos($t) + 1);
my $D64p3 =  exp( 4*$i*$f) * sqrt($X*$F)/32   * sin($t)    * (1+cos($t))**3 * (11*cos($t)**2 - 11*cos($t) + 2);
my $D64m3 =  exp(-4*$i*$f) * sqrt($X*$F)/32   * sin($t)    * (1-cos($t))**3 * (11*cos($t)**2 + 11*cos($t) + 2);
my $D64p4 =  exp( 4*$i*$f) * (1/32)           *              (1+cos($t))**4 * (33*cos($t)**2 - 44*cos($t) + 13);
my $D64m4 =  exp(-4*$i*$f) * (1/32)           *              (1-cos($t))**4 * (33*cos($t)**2 + 44*cos($t) + 13);
my $D64p5 =  exp( 4*$i*$f) * sqrt($T*$E)/32   * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D64m5 = -exp(-4*$i*$f) * sqrt($T*$E)/32   * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D64p6 =  exp( 4*$i*$f) * sqrt($X*$E)/64   * sin($t)**2 * (1+cos($t))**4;
my $D64m6 =  exp(-4*$i*$f) * sqrt($X*$E)/64   * sin($t)**2 * (1-cos($t))**4;

my $D65p0 =  exp( 5*$i*$f) * 3*sqrt($S*$E)/16  * sin($t)**5 * cos($t);
my $D65p1 =  exp( 5*$i*$f) * sqrt($X*$E)/32    * sin($t)**4 * (1+cos($t))    * (6*cos($t) - 1);
my $D65m1 =  exp(-5*$i*$f) * sqrt($X*$E)/32    * sin($t)**4 * (1-cos($t))    * (6*cos($t) + 1);
my $D65p2 =  exp( 5*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**3 * (1+cos($t))**2 * (3*cos($t) - 1);
my $D65m2 =  exp(-5*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**3 * (1-cos($t))**2 * (3*cos($t) + 1);
my $D65p3 =  exp( 5*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**2 * (1+cos($t))**3 * (2*cos($t) - 1);
my $D65m3 =  exp(-5*$i*$f) * sqrt($e*$F*$E)/32 * sin($t)**2 * (1-cos($t))**3 * (2*cos($t) + 1);
my $D65p4 =  exp( 5*$i*$f) * sqrt($T*$E)/32    * sin($t)    * (1+cos($t))**4 * (3*cos($t) - 2);
my $D65m4 =  exp(-5*$i*$f) * sqrt($T*$E)/32    * sin($t)    * (1-cos($t))**4 * (3*cos($t) + 2);
my $D65p5 =  exp( 5*$i*$f) * (1/32)            *              (1+cos($t))**5 * (6*cos($t) - 5);
my $D65m5 =  exp(-5*$i*$f) * (1/32)            *              (1-cos($t))**5 * (6*cos($t) + 5);
my $D65p6 =  exp( 5*$i*$f) * sqrt($e)/32       * sin($t)    * (1+cos($t))**5;
my $D65m6 = -exp(-5*$i*$f) * sqrt($e)/32       * sin($t)    * (1-cos($t))**5;

my $D66p0 =  exp( 6*$i*$f) * sqrt($e*$S*$E)/32 * sin($t)**6;
my $D66p1 =  exp( 6*$i*$f) * 3*sqrt($T*$E)/32  * sin($t)**5 * (1+cos($t));
my $D66m1 =  exp(-6*$i*$f) * 3*sqrt($T*$E)/32  * sin($t)**5 * (1-cos($t));
my $D66p2 =  exp( 6*$i*$f) * 3*sqrt($F*$E)/64  * sin($t)**4 * (1+cos($t))**2;
my $D66m2 =  exp(-6*$i*$f) * 3*sqrt($F*$E)/64  * sin($t)**4 * (1-cos($t))**2;
my $D66p3 =  exp( 6*$i*$f) * sqrt($F*$E)/32    * sin($t)**3 * (1+cos($t))**3;
my $D66m3 =  exp(-6*$i*$f) * sqrt($F*$E)/32    * sin($t)**3 * (1-cos($t))**3;
my $D66p4 =  exp( 6*$i*$f) * sqrt($X*$E)/64    * sin($t)**2 * (1+cos($t))**4;
my $D66m4 =  exp(-6*$i*$f) * sqrt($X*$E)/64    * sin($t)**2 * (1-cos($t))**4;
my $D66p5 =  exp( 6*$i*$f) * sqrt($e)/32       * sin($t)    * (1+cos($t))**5;
my $D66m5 =  exp(-6*$i*$f) * sqrt($e)/32       * sin($t)    * (1-cos($t))**5;
my $D66p6 =  exp( 6*$i*$f) * (1/64)            *              (1+cos($t))**6;
my $D66m6 =  exp(-6*$i*$f) * (1/64)            *              (1-cos($t))**6;

my $D6 = [
     [ $D66p6,  $D66p5,  $D66p4,  $D66p3,  $D66p2,  $D66p1,  $D66p0,  $D66m1,  $D66m2,  $D66m3,  $D66m4,  $D66m5,  $D66m6],
     [-$D65p6,  $D65p5,  $D65p4,  $D65p3,  $D65p2,  $D65p1,  $D65p0,  $D65m1,  $D65m2,  $D65m3,  $D65m4,  $D65m5, -$D65m6],
     [ $D64p6, -$D64p5,  $D64p4,  $D64p3,  $D64p2,  $D64p1,  $D64p0,  $D64m1,  $D64m2,  $D64m3,  $D64m4, -$D64m5,  $D64m6],
     [-$D63p6,  $D63p5, -$D63p4,  $D63p3,  $D63p2,  $D63p1,  $D63p0,  $D63m1,  $D63m2,  $D63m3, -$D63m4,  $D63m5, -$D63m6],
     [ $D62p6, -$D62p5,  $D62p4, -$D62p3,  $D62p2,  $D62p1,  $D62p0,  $D62m1,  $D62m2, -$D62m3,  $D62m4, -$D62m5,  $D62m6],
     [-$D61p6,  $D61p5, -$D61p4,  $D61p3, -$D61p2,  $D61p1,  $D61p0,  $D61m1, -$D61m2,  $D61m3, -$D61m4,  $D61m5, -$D61m6],
     [ $D60p6, -$D60p5,  $D60p4, -$D60p3,  $D60p2, -$D60p1,  $D60p0,  $D60p1,  $D60p2,  $D60p3,  $D60p4,  $D60p5,  $D60p6],
     [ $D61m6,  $D61m5,  $D61m4,  $D61m3,  $D61m2,  $D61m1, -$D61p0,  $D61p1,  $D61p2,  $D61p3,  $D61p4,  $D61p5,  $D61p6],
     [ $D62m6,  $D62m5,  $D62m4,  $D62m3,  $D62m2, -$D62m1,  $D62p0, -$D62p1,  $D62p2,  $D62p3,  $D62p4,  $D62p5,  $D62p6],
     [ $D63m6,  $D63m5,  $D63m4,  $D63m3, -$D63m2,  $D63m1, -$D63p0,  $D63p1, -$D63p2,  $D63p3,  $D63p4,  $D63p5,  $D63p6],
     [ $D64m6,  $D64m5,  $D64m4, -$D64m3,  $D64m2, -$D64m1,  $D64p0, -$D64p1,  $D64p2, -$D64p3,  $D64p4,  $D64p5,  $D64p6],
     [ $D65m6,  $D65m5, -$D65m4,  $D65m3, -$D65m2,  $D65m1, -$D65p0,  $D65p1, -$D65p2,  $D65p3, -$D65p4,  $D65p5,  $D65p6],
     [ $D66m6, -$D66m5,  $D66m4, -$D66m3,  $D66m2, -$D66m1,  $D66p0, -$D66p1,  $D66p2, -$D66p3,  $D66p4, -$D66p5,  $D66p6]
          ];

my ($O6m6, $O6m5, $O6m4, $O6m3, $O6m2, $O6m1) = symbols(qw(Osms OsmF Osmf OsmT Osmt Osmo));
my $O6p0 = symbols(qw(Ospz));
my ($O6p6, $O6p5, $O6p4, $O6p3, $O6p2, $O6p1) = symbols(qw(Osps OspF Ospf OspT Ospt Ospo));


      
# ------------------------------------ Finished !! ------------------------------- #

# The operator in the rotated frame. Rx(y) is equivalent to {Oxy} in Rudowicz notation
# and Ox_y, Oxmy is equivalent to [Oxy], [OxyM] in Rudowicz

my $B6 = [ [$O6m6],[$O6m5],[$O6m4],[$O6m3],[$O6m2],[$O6m1],[$O6p0],[$O6p1],[$O6p2],[$O6p3],[$O6p4],[$O6p5],[$O6p6] ];

my $tmp1 = mmult($invA6,$B6);
my $tmp2 = mmult($D6, $tmp1);
my $R6 = mmult($A6, $tmp2);

my $S6 = [];
my $M6 = [];
for $l (0 .. 12) { 
  $S6->[$l] = $R6->[$l][0]->sub(f=>$phi_val, t=>$theta_val);
  $S6->[$l] = $S6->[$l]->sub(two=>2,three=>3,five=>5,six=>6,seven=>7,eleven=>11); 
  $M6->[0][$l] = $S6->[$l]->sub(Osms=>1,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[1][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>1,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[2][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>1,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[3][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>1,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[4][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>1,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[5][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>1,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[6][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>1,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[7][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>1,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[8][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>1,OspT=>0,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[9][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>1,Ospf=>0,OspF=>0,Osps=>0);
  $M6->[10][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>1,OspF=>0,Osps=>0);
  $M6->[11][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>1,Osps=>0);
  $M6->[12][$l] = $S6->[$l]->sub(Osms=>0,OsmF=>0,Osmf=>0,OsmT=>0,Osmt=>0,Osmo=>0,Ospz=>0,Ospo=>0,Ospt=>0,OspT=>0,Ospf=>0,OspF=>0,Osps=>1);
}

if (!$output) {
  %T6 = (0,'O6m6',1,'O6m5',2,'O6m4',3,'O6m3',4,'O6m2',5,'O6m1',6,'O6_0',7,'O6_1',8,'O6_2',9,'O6_3',10,'O6_4',11,'O6_5',12,'O6_6');
  print "Crystal field parameters under rotation of phi=$phi about y-axis and theta=$theta about z-axis are given by:\n\n";
  print "Transformations relations for O6q:\n";
  for $l (0 .. 12) { 
    print "\[$T6{$l}\] = "; 
    $_ = $S6->[$l]; 
    s/\$Osms/\{O6m6\}/g;
    s/\$OsmF/\{O6m5\}/g;
    s/\$Osmf/\{O6m4\}/g;
    s/\$OsmT/\{O6m3\}/g;
    s/\$Osmt/\{O6m2\}/g;
    s/\$Osmo/\{O6m1\}/g;
    s/\$Ospz/\{O6_0\}/g;
    s/\$Ospo/\{O6_1\}/g;
    s/\$Ospt/\{O6_2\}/g;
    s/\$OspT/\{O6_3\}/g;
    s/\$Ospf/\{O6_4\}/g;
    s/\$OspF/\{O6_5\}/g;
    s/\$Osps/\{O6_6\}/g;
    print "$_\n"; 
  }
  print "\nRotation Matrix for O6q:\n";
  for $l (0 .. 12) { for $m (0 .. 12) {
    print $M6->[$l][$m]; print "\t"; } print "\n"; }
}
else {
  %T6 = (0,'O_6^{-6}',1,'O_6^{-5}',2,'O_6^{-4}',3,'O_6^{-4}',4,'O_6^{-2}',5,'O_6^{-1}',6,'O_6^0',7,'O_6^1',8,'O_6^2',9,'O_6^3',10,'O_6^4',11,'O_6^5',12,'O_6^6');

  if (!$z_int) {
    open (outfile, ">$output") or die "$0: cannot open $output for output.";
    print outfile << "EOF";

\\documentclass [12pt,a4paper,notitlepage]{article}
\\usepackage{latexsym}
\\usepackage{pslatex}
\\setlength{\\unitlength}{1cm}

\\begin{document}

Crystal field parameters under rotation of phi=$phi\$^{o}\$ about y-axis and theta=$theta\$^{o}\$ about z-axis are given by:

Transformation relations for \$O_6^q\$ Steven's equivalent operators:

EOF
  
  }
  else {
    open (outfile, ">>$output") or die;
    print outfile "\nTransformation relations for \$O_6^q\$ Steven's equivalent operators:\n";
  }

  for $l (0 .. 12) {
    $_ = $S6->[$l];
    s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
    s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;
    s/\$Osms/\\left\[ O\_6\^\{-6\} \\right\]/g;
    s/\$OsmF/\\left\[ O\_6\^\{-5\} \\right\]/g;
    s/\$Osmf/\\left\[ O\_6\^\{-4\} \\right\]/g;
    s/\$OsmT/\\left\[ O\_6\^\{-3\} \\right\]/g;
    s/\$Osmt/\\left\[ O\_6\^\{-2\} \\right\]/g;
    s/\$Osmo/\\left\[ O\_6\^\{-1\} \\right\]/g;
    s/\$Ospz/\\left\[ O\_6\^0 \\right\]/g;
    s/\$Ospo/\\left\[ O\_6\^1 \\right\]/g;
    s/\$Ospt/\\left\[ O\_6\^2 \\right\]/g;
    s/\$OspT/\\left\[ O\_6\^3 \\right\]/g;
    s/\$Ospf/\\left\[ O\_6\^4 \\right\]/g;
    s/\$OspF/\\left\[ O\_6\^5 \\right\]/g;
    s/\$Osps/\\left\[ O\_6\^6 \\right\]/g;
    s/\*//g;
    print outfile "\\begin\{equation\}\n";
    print outfile "\\left\\\{ $T6{$l} \\right\\\} = $_\n";
    print outfile "\\end\{equation\}\n";
  }

  print outfile "\nRotation Matrix for \$O_6^q\$:\n\\begin{equation}\n\\left(\n\\begin{array}{ccccccccccccc}\n";
  for $l (0 .. 12) { 
    for $m (0 .. 12) {
      $_ = $M6->[$l][$m];
      s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
      s/\/(\\sqrt\{[\/\.\d]*\})/\\frac\{\}\{$1\}/g;
      s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;       # + matches 1 or more times
      s/\*//g;
      if ($m == 12) { print outfile "$_ "; } else {
      print outfile "$_ & "; }
    } 
    print outfile "\\\\\n"; 
  }
  print outfile "\\end{array}\n\\right)\n\\end{equation}\n";

  print outfile "\\end\{document\}\n";

  close (out_file);
}