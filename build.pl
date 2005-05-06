#!/home/ivan/bin/perl -s

use strict;
use warnings;

use lib 'lib';
use Chemistry::3DBuilder qw(:all);
use Chemistry::File::SMILES;
use Chemistry::File::MDLMol;
use Chemistry::File::Mopac;

our $DEBUG;
our $format;
my $smiles = shift || die "pls give a SMILES string\n";
my $mol = Chemistry::Mol->parse($smiles, format => 'smiles');

$Chemistry::3DBuilder::DEBUG = $DEBUG;
build_3d($mol);

print $mol->print(format => $format || 'mdl');
