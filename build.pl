#!/home/ivan/bin/perl -s

use strict;
use warnings;

use lib 'lib';
use Chemistry::3DBuilder;
use Chemistry::File::SMILES;
use Chemistry::File::MDLMol;

our $DEBUG;
my $mol = Chemistry::Mol->parse(shift || 'CCC', format => 'smiles');

$Chemistry::3DBuilder::DEBUG = $DEBUG;
Chemistry::3DBuilder::build3d($mol);

print $mol->print(format => 'mdl');
