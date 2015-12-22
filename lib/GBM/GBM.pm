package GBM::GBM;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Dancer2;


our $VERSION = '0.1';
my $registry = 'Bio::EnsEMBL::Registry';

get '/' => sub {
    
    $registry->load_all('./loutre_config');
   
print Dumper $registry; 
    template 'index', {};
};

true;
