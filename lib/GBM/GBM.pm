package GBM::GBM;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Cwd;
use File::Slurp;
use Dancer2;


our $VERSION = '0.1';
my $registry = 'Bio::EnsEMBL::Registry';
my $cw_dir = getcwd;
$registry->load_all("$cw_dir/environments/loutre_config");

get '/' => sub {
  
 my @species = @{ $registry->get_all_species };
 my (%trans_biotypes, %sources);
 load_data(\%trans_biotypes, \@species, 'biotypes');
 load_data(\%sources,  \@species, 'sources');
# foreach my $species( @species ) {
#  my $trans_adp = $registry->get_adaptor( "$species", 'core', 'transcript' );
#  foreach my $transcript( @{ $trans_adp->fetch_all } ) {
#   $trans_biotypes{ $transcript->biotype }{ $species }++;
#  }
# } 
 
print Dumper \%trans_biotypes;
 
 template 'index', {
  'species_list' => \@species,
  'biotypes'     => \%trans_biotypes,
  'sources'      => \%sources,

 };

};

sub load_data {
 my ($data_hash, $species_list, $data_type) = @_;
 foreach my $species(@{ $species_list }) {
  my $file_path = "$cw_dir/data/$data_type.$species";
  if( -e $file_path ) { 
   foreach my $data_type(read_file($file_path)) {
    chomp $data_type;
    push @{ $data_hash->{ $data_type } }, $species;
   }
  }
 }
}

true;
