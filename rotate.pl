#!/usr/bin/env perl

use strict;
use feature 'say';
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Getopt::ArgParse;
use Scalar::Util 'looks_like_number';
#https://metacpan.org/pod/pp

=sub dihedral_angle {
	my ($c1, $c2, $c3, $c4) = @_; # c means "coordinate"

# Calculate the vectors between the points
	my @b1 = (
					$c2->[0] - $c1->[0],	# x coordinate
					$c2->[1] - $c1->[1], # y coordinate
					$c2->[2] - $c1->[2]  # z coordinate
				);
	
	my @b2 = (
					$c3->[0] - $c2->[0],	# x coordinate
					$c3->[1] - $c2->[1], # y coordinate
					$c3->[2] - $c2->[2]  # z coordinate
				);
	my @b3 = (
					$c4->[0] - $c3->[0],	# x coordinate
					$c4->[1] - $c3->[1], # y coordinate
					$c4->[2] - $c3->[2]  # z coordinate
	);
	my $magnitude = sqrt($b1[0]**2 + $b1[1]**2 + $b1[2]**2);
	@b1 = map {$_ / $magnitude} @b1;
	$magnitude = sqrt($b2[0]**2 + $b2[1]**2 + $b2[2]**2);
	@b2 = map {$_ / $magnitude} @b2;
	$magnitude = sqrt($b3[0]**2 + $b3[1]**2 + $b3[2]**2);
	@b3 = map {$_ / $magnitude} @b3;
	$magnitude = 0;
#step 2: Cross Products
   my @n1 = ($b1[1]*$b2[2]-$b1[2]*$b2[1],$b1[2]*$b2[0]-$b1[0]*$b2[2],$b1[0]*$b2[1]-$b1[1]*$b2[0]);
   for (my $i = 0; $i < @n1; $i++) {
      $magnitude += ($n1[$i])**2;
   }
   $magnitude = $magnitude**.5;
   for (my $i = 0; $i < @n1; $i++) {
      $n1[$i] = $n1[$i]/$magnitude;
   }
   $magnitude = 0;
   my @n2 = ($b2[1]*$b3[2]-$b2[2]*$b3[1],$b2[2]*$b3[0]-$b2[0]*$b3[2],$b2[0]*$b3[1]-$b2[1]*$b3[0]);
   for (my $i = 0; $i < @n2; $i++) {
      $magnitude += ($n2[$i])**2;
   }
   $magnitude = $magnitude**.5;
   for (my $i = 0; $i < @n2; $i++) {
      $n2[$i] = $n2[$i]/$magnitude;
   }
   $magnitude = 0;
#step 3: make x and y; m1 = n1 x b2
   my @m1 = ($n1[1]*$b2[2]-$n1[2]*$b2[1],$n1[2]*$b2[0]-$n1[0]*$b2[2],$n1[0]*$b2[1]-$n1[1]*$b2[0]);
   my $x = $n1[0]*$n2[0]+$n1[1]*$n2[1]+$n1[2]*$n2[2];
   my $y = $m1[0]*$n2[0]+$m1[1]*$n2[1]+$m1[2]*$n2[2];
#step 4: atan2
   return (-180/atan2(0,-1))*atan2($y,$x);
}
=cut

sub rotate_point {
	my ($p_ref, $a_ref, $b_ref, $theta) = @_;
	$theta *= atan2(0,-1)/180; # convert to radians
	# Calculate the vector components for the axis of rotation
	my $ux = $b_ref->[0] - $a_ref->[0];
	my $uy = $b_ref->[1] - $a_ref->[1];
	my $uz = $b_ref->[2] - $a_ref->[2];
	my $u_norm = sqrt($ux**2 + $uy**2 + $uz**2);
	$ux /= $u_norm;
	$uy /= $u_norm;
	$uz /= $u_norm;

	# Calculate the sine and cosine of the rotation angle
	my $cos_theta = cos($theta);
	my $sin_theta = sin($theta);

	# Calculate the dot product of P and the axis of rotation
	my $px = $p_ref->[0];
	my $py = $p_ref->[1];
	my $pz = $p_ref->[2];
	my $dot_product = $px*$ux + $py*$uy + $pz*$uz;

	# Calculate the vector components for the parallel and perpendicular components of P
	my $parallel_x = $ux*$dot_product;
	my $parallel_y = $uy*$dot_product;
	my $parallel_z = $uz*$dot_product;
	my $perpendicular_x = $px - $parallel_x;
	my $perpendicular_y = $py - $parallel_y;
	my $perpendicular_z = $pz - $parallel_z;

	# Calculate the rotated vector components
	my $rotated_x = $parallel_x + ($perpendicular_x*$cos_theta) + (($uy*$perpendicular_z) - ($uz*$perpendicular_y))*$sin_theta;
	my $rotated_y = $parallel_y + ($perpendicular_y*$cos_theta) + (($uz*$perpendicular_x) - ($ux*$perpendicular_z))*$sin_theta;
	my $rotated_z = $parallel_z + ($perpendicular_z*$cos_theta) + (($ux*$perpendicular_y) - ($uy*$perpendicular_x))*$sin_theta;

	return [$rotated_x, $rotated_y, $rotated_z];
}

my $parser = Getopt::ArgParse->new_parser(
	help   => 'This script rotates atoms in an XYZ coordinate file.  All atoms are 0-indexed.',
	epilog => 'Example: perl rotate.pl -i ethane.xyz -a1 2 -a2 0 -o ethane.rotated.x.xyz -rotated 1,3,4 -th 60'
);
$parser-> add_args(
	['-atom1',   '-a1', required => 1, type => 'Scalar', help => 'the 1st atom defining the axis of rotation'],
	['-atom2',   '-a2', required => 1, type => 'Scalar', help => 'the 2nd atom defining the axis of rotation'],
	['-input',   '-i',  required => 1, type => 'Scalar', help => 'the source file to be read'],
	['-output',  '-o',  required => 1, type => 'Scalar', help => 'output file'],
	['-rotated', '-r',  required => 1, type => 'Scalar', help => 'a comma-delimited list of atoms that need to be rotated'],
	['-theta',   '-t',  required => 1, type => 'Scalar', help => 'the angle that the atom should be rotated']
);

my $ns = $parser->parse_args(@ARGV);
# check inputs
if ($ns->atom1 !~ m/^\d+$/) {
	my $a1 = $ns->atom1;
	die "atom1 must be in [0,inf] but you entered \"$a1\""
}
if ($ns->atom2 !~ m/^\d+$/) {
	my $a2 = $ns->atom2;
	die "atom2 must be in [0,inf] but you entered \"$a2\""
}
unless (looks_like_number($ns->theta)) {
	my $th = $ns->theta;
	die "theta must be numeric, but you entered $th"
}
my @rotated_atoms = split ',', $ns->rotated;
my @non_numeric_rotate = grep {$_ !~ m/^(\d+)$/} @rotated_atoms;

if (scalar @non_numeric_rotate > 0) {
	p @non_numeric_rotate;
	die '"rotated atoms" must be a list of numbers in [0,inf], the above values cannot work'
}
if (grep {$_ == $ns->atom1 || $_ == $ns->atom2} @rotated_atoms) {
	die "Neither a1 nor a2 can be included in rotated atoms.";
}
open my $in,  '<', $ns->input;
open my $out, '>', $ns->output;
my (@elements, @coords);
while (<$in>) {
	if ($. <= 2) { # I'm not error checking for the first 2 lines
		print $out $_;
		next
	}
	my @line = split;
	if (scalar @line != 4) {
		p @line;
		die 'the current line may have a formatting error, there should be exactly 4 elements.';
	}
	my @non_numeric = grep {not looks_like_number($_)} @line[1..3];
	if (scalar @non_numeric > 0) {
		print STDERR "from line $.\n";
		p @line;
		print STDERR "the below values are non-numeric:\n";
		p @non_numeric;
		die 'all values must be numeric'
	}
	push @elements, shift @line;
	push @coords,   [@line];
}
close $in;

my ($a1, $a2) = ($ns->atom1, $ns->atom2);
my $n = scalar @coords - 1;
if ($ns->atom1 > scalar @coords) {
	die "a1 = $a1 doesn't exist, as the max index is $n"
}
if ($ns->atom2 > scalar @coords) {
	die "a2 = $a2 doesn't exist, as the max index is $n"
}
my @atoms_above_limit = grep { $_ > $n } @rotated_atoms;
if (scalar @atoms_above_limit > 0) {
	p @atoms_above_limit;
	die "the above atom numbers are not available, as the max atom # is $n"
}
#say dihedral_angle($coords[1], $coords[0], $coords[2], $coords[5]);
foreach my $ra (@rotated_atoms) {
#	my $dihedral1 = calculate_dihedral(@{ $d);
	$coords[$ra] = rotate_point($coords[$ra], $coords[$ns->atom1], $coords[$ns->atom2], $ns->theta);
}
#say dihedral_angle($coords[1], $coords[0], $coords[2], $coords[5]);
foreach my $n (0..$#coords) {
	print $out join (' ', $elements[$n], @{ $coords[$n] }) . "\n";
}
close $out;
my $outfile = $ns->output;
print "Just wrote $outfile\n";
