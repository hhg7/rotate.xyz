# rotate.xyz
Rotate a group of atoms in a XYZ coordinate file for computational chemistry
rotate.pl: This script rotates atoms in an XYZ coordinate file.  All atoms
are 0-indexed.

usage: rotate.pl -theta|-t -rotated|-r -output|-o -input|-i -atom2|-a2
-atom1|-a1

required named arguments:
  -theta, -t THETA          the angle that the atom should be rotated
  -rotated, -r ROTATED      a comma-delimited list of atoms that need to be rotated
  -output, -o OUTPUT        output file
  -input, -i INPUT          the source file to be read
  -atom2, -a2 ATOM2         the 2nd atom defining the axis of rotation
  -atom1, -a1 ATOM1         the 1st atom defining the axis of rotation

Example: perl rotate.pl -i ethane.xyz -a1 2 -a2 0 -o ethane.rotated.x.xyz
-rotated 1,3,4 -th 60

requires Data::Printer, Getopt::ArgParse, and Devel::Confess
