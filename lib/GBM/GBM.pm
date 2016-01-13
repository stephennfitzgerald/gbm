package GBM::GBM;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use String::Random qw(random_regex);
use Cwd;
use Carp;
use File::Slurp;
use Dancer2;

our $VERSION = '0.1';
our $download_dir = './public/downloads/';
our $download_file_path = 'downloads/';
my $cw_dir = getcwd;

our %assembly_hack = ( # a hack in order to get the correct seq_regions
 human => '-38',
 mouse => '-38',
); 

our %file_types = (
 bed    => \&bed,
# gtf    => \&gtf,
 fasta  => \&fasta,
);

get '/' => sub {
  
 my $registry = 'Bio::EnsEMBL::Registry';
 $registry->load_all("$cw_dir/environments/loutre_config");
 my @species_list = @{ $registry->get_all_species };
 my (%sources, %trans_biotypes, @transcripts, $no_results, $download_file, $slice_adaptor, $file_ok);
 load_data(\%trans_biotypes, \@species_list, 'biotypes');
 load_data(\%sources,  \@species_list, 'sources');

 my ($slices, $seq_reg, $chrom, $start, $end);
 my $species = param('species');
 my $genomic_region = param('genomic_region');
 my $biotype = param('feature_biotype');
 my $source  = param('feature_source');
 my $strand = param('genomic_strand');
 if($strand && $strand eq 'both') {
  $strand = 0;
 }
 elsif($strand) {
  $strand = $strand eq 'forward' ? 1 : -1;
 }

 if($genomic_region && $species) {
  $slice_adaptor = $registry->get_adaptor( "$species", "core", "slice" );
  ($seq_reg, $start, $end) = $genomic_region=~/\s*(\w+)\s*[:-]+\s*(\d+)\s*-+\s*(\d+)\s*$/;
  $seq_reg=~s/^chr//;
  $seq_reg .= $assembly_hack{ "$species" };
  $chrom = 'chr' . $seq_reg;
  my $slice = $slice_adaptor->fetch_by_region('toplevel', "$seq_reg", $start, $end);
  if($slice) {
   $slices = [ $slice ];
  }
  else {
   $slices = [ $slice_adaptor->fetch_by_region('toplevel', "$chrom", $start, $end) ];
  } 
 }
 elsif($species) {
  $slice_adaptor = $registry->get_adaptor( "$species", "core", "slice" );
  $slices = $slice_adaptor->fetch_all('toplevel');
 }
 if(! $slices->[0] && $species) {
  croak "Can't find genomic regions $chrom:$start-$end for species $species :";
 } 

 foreach my $slice( @$slices ) {
  last if (($source eq 'all' && param('not_source')) || ($biotype eq 'all' && param('not_biotype')));
  foreach my $transcript(@{ $slice->get_all_Transcripts }) {
   next if (! $transcript->is_current);
   if($strand) {
    next if($transcript->strand ne $strand);
   }
   if($biotype ne 'all') {
    if(param('not_biotype')) {
     next if($transcript->biotype eq $biotype);
    }
    else {
      next if($transcript->biotype ne $biotype);
    }
   }
   if($source ne 'all') {
    if(param('not_source')) {
     next if($transcript->source eq $source);
    }
    else {
     next if($transcript->source ne $source);
    }
   }
   push@transcripts, $transcript;
  }
 } 

 if(! scalar @transcripts && $species) {
  $no_results = 1;
 }
 elsif(scalar @transcripts) {
  my $file_type = param('feature_format');
  my $file_name = 'feat_' . random_regex('\w'x15) . ".$file_type.txt";
  my $file_loc  = $download_dir . $file_name;
  open my $fh, '>', $file_loc or croak "Can't open file $file_name\n";
  if($file_types{$file_type}->($fh, \@transcripts)) {
   $file_ok = $download_file_path . $file_name;
  }
 }

 template 'index', {
  'species_list'  => \@species_list,
  'biotypes'      => \%trans_biotypes,
  'sources'       => \%sources,
  'no_results'    => $no_results,
  'download_file' => $download_file,
  'file_types'    => \%file_types,
  'file_loc'      => $file_ok,
 };

};

sub bed {
 my($fh, $transcripts) = @_;
 foreach my $trans(@{ $transcripts }) {
  my ($seq_region_name) = $trans->seq_region_name=~/(\w+)-\d+$/xms;
  print $fh join("\t", $seq_region_name, ($trans->seq_region_start - 1), $trans->seq_region_end, $trans->seq_region_strand, 
            $trans->stable_id, $trans->biotype, $trans->source, $trans->get_Gene()->stable_id), "\n";
 }
 close $fh;
 return 1;
}

#sub gtf {
# my($fh, $transcripts) = @_; 
# foreach my $trans(@{ $transcripts }) {
# }
#}

sub fasta {
 my($fh, $transcripts) = @_;
 foreach my $trans(@{ $transcripts }) {
  my ($seq_region_name) = $trans->seq_region_name=~/(\w+)-\d+$/xms;
  print $fh '>', $seq_region_name, q{:}, $trans->seq_region_start, q{-}, $trans->seq_region_end, q{:}, $trans->seq_region_strand,
            q{::}, join( q{::}, $trans->stable_id, $trans->biotype, $trans->source, $trans->get_Gene()->stable_id), "\n",
            $trans->seq()->seq, "\n";
 }
 close $fh;
 return 1;
}

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
