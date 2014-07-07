#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $version = "0.0.23";

# version update
# 0.0.23
#   -added option to require a minimum level of evidence for sample level genotyping
#
# 0.0.22
#   -updated version to be consistent with dinumt package
#
# 0.0.16
#   -fixed program with too many supporting reads crashing GL calculation
#   -updated output to include more original INFO tag information
#
# 0.0.15
#   -added option to turn off EM iterations
#
# 0.0.14
#   -added MT mapping check for clipped sequences at potential breakpoints
#
# 0.0.13
#   -fixed error with likelihood calculation
#   -added support for by_chr_dir genotyping
#
# 0.0.12
#   -added per-read mapQ errors with likelihood calculation
#
# 0.0.11
#   -updated supporting read determination
#
# 0.0.10
#   -integrated RP and SR liklihoods, when available
#
# 0.0.9
#   -incorporated EM algorithm for implementing allele frequency priors
#
# 0.0.8
#   -updated support counting around breakpoints
#
# 0.0.7
#   -modified scoring to implement bayesian scoring of evidence
#   -incorporated insert length cutoff for alternative allele determination
#
# 0.0.6
#   -fixed bug where altSR was not being initialized correctly
#   -removed max_read_cov filter in report()
#
# 0.0.5
#   -added --chr option for by-chromosome analysis
#
# 0.0.4
#   -bug fixes
#
# 0.0.3
#   -modified scoring function to better assess clipping in breakpoint determination and scoring
#   -modified rp scoring to not consider reads clipped at breakpoint as reference
#
# 0.0.2
#   -modified input parameters to expect VCF format
#
# 0.0.1
#   -initial version

my %opts = ();

$opts{len_cluster_include} = 600;
$opts{len_cluster_link}    = 800;
$opts{min_reads_cluster}   = 1;
$opts{min_clipped_seq}     = 5;
$opts{clipped_flank}       = 50;
$opts{max_num_clipped}     = 5;
$opts{include_mask}        = 0;
$opts{min_evidence}        = 3;
$opts{min_map_qual}        = 10;
$opts{filter_qual}         = 13;
$opts{filter_depth}        = 5;
$opts{max_read_cov}        = 200;
$opts{use_priors}          = 0;
$opts{mask_filename}       = "/home2/remills/projects/numts/numtS.bed";
$opts{info_filename}       = "/scratch/remills_flux/remills/sampleInfo.txt";
$opts{reference}           = "/nfs/remills-scratch/reference/hs37d5/hs37d5.fa";
$opts{mt_filename}         = "/nfs/remills-scratch/reference/GRCh37/Sequence/Chromosomes/MT.fa";
$opts{samtools}            = "/home2/remills/bin/samtools";
$opts{exonerate}           = "/home2/remills/apps/exonerate-2.2.0/bin/exonerate";
$opts{dir_tmp}             = "/tmp";
$opts{prefix}              = "numt";
$opts{len_mt}              = 16596;                                                                #eventually should be read in by BAM header
$opts{ploidy}              = 2;

my $optResult = GetOptions(
    "input_filename=s"      => \$opts{input_filename},
    "output_filename=s"     => \$opts{output_filename},
    "mask_filename=s"       => \$opts{mask_filename},
    "info_filename=s"       => \$opts{info_filename},
    "mt_filename=s"         => \$opts{mt_filename},
    "dir_tmp=s"             => \$opts{dir_tmp},
    "chr=s"                 => \$opts{chr},
    "include_mask"          => \$opts{include_mask},
    "len_cluster_include=i" => \$opts{len_cluster_include},
    "len_cluster_link=i"    => \$opts{len_cluster_link},
    "min_reads_cluster=i"   => \$opts{min_reads_cluster},
    "min_evidence=i"        => \$opts{min_evidence},
    "min_clipped_seq=i"     => \$opts{min_clipped_seq},
    "max_num_clipped=i"     => \$opts{max_num_clipped},
    "min_map_qual=i"        => \$opts{min_map_qual},
    "max_read_cov=i"        => \$opts{max_read_cov},
    "read_groups"           => \$opts{read_groups},
    "breakpoint"            => \$opts{breakpoint},
    "reference=s"           => \$opts{reference},
    "samtools=s"            => \$opts{samtools},
    "exonerate=s"           => \$opts{exonerate},
    "use_priors"            => \$opts{use_priors},
    "by_chr_dir"            => \$opts{by_chr_dir},
    "prefix=s"              => \$opts{prefix},
    "ucsc"                  => \$opts{ucsc},
    "help"                  => \$opts{help},
    "verbose"               => \$opts{verbose}
);

#checkOptions( $optResult, \%opts, $version );

my $seq_num  = 0;
my %seq_hash = ();

my %sorted_hash = ();

my $i            = 1;
my %sample_hash  = ();
my %infile_hash  = ();
my %group_hash   = ();
my %outfile_hash = ();
my %mask_hash    = ();
my %data_hash    = ();

getInput( \%infile_hash, \%mask_hash, \%sample_hash );
getData( \%infile_hash, \%mask_hash, \%sample_hash, \%data_hash );
assessBreaks( \%infile_hash, \%sample_hash, \%data_hash );
refineData( \%infile_hash, \%mask_hash, \%sample_hash, \%data_hash );
scoreData( \%data_hash, \%sample_hash );
report( \%infile_hash, \%data_hash );

################################################################################################################
sub mappedClippedSeq {
    my ($clipSeq) = @_;

    if ( $clipSeq =~ /N/ ) { return 0; }

    my $random    = int( rand(100000) );
    my @timeData  = localtime(time);
    my $timestamp = join( '', @timeData );
    my $fileClp   = "$opts{dir_tmp}/clp$timestamp$random.fa";
    open( CLP, ">$fileClp" );
    print CLP ">clp\n$clipSeq\n";
    close CLP;

    my $results = `$opts{exonerate} --model affine:local --exhaustive yes --percent 70 --showvulgar no --showalignment no --showcigar yes $fileClp $opts{mt_filename} 2> /dev/null`;
    my ( $query_id, $query_start, $query_end, $query_strand, $target_id, $target_start, $target_end, $target_strand, $score, $cigar ) = $results =~ /cigar\:\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\+|-)\s+(\d+)\s+(.*?)\n/;
    `rm $fileClp`;
    if ($cigar) { return 1; }
    else { return 0; } 
}

sub assessBreaks {
    my ( $infile_hash, $sample_hash, $data_hash ) = @_;
    print "Entering assessBreaks()\n" if $opts{verbose};
    foreach my $var ( keys %$data_hash ) {
        print "-variant $var\n" if $opts{verbose};
        my %clippedPos = ();
        my %clippedSeq = ();
        my $sumClipped = 0;
        my $winStart   = 1e10;
        my $winEnd     = 0;
        my $maxFirst   = 0;
        my $maxSecond  = 0;

        foreach my $sample ( keys %{ $data_hash{$var} } ) {
            foreach my $cPos ( keys %{ $data_hash{$var}{$sample}{clipPos} } ) {
                next if $cPos == -1;
                if ( defined( $data_hash{$var}{$sample}{clipPos}{$cPos}{1} ) ) {
                    if ( mappedClippedSeq( $data_hash{$var}{$sample}{clipSeq}{$cPos} ) ) {
                        $clippedPos{$cPos} += $data_hash{$var}{$sample}{clipPos}{$cPos}{1};
                        $sumClipped += $data_hash{$var}{$sample}{clipPos}{$cPos}{1};

                        #if (!defined($clippedSeq{$cPos} || length($data_hash{$var}{$sample}{clipSeq}{$cPos}) > $clippedSeq{$cPos})) { $clippedSeq{$cPos} = $data_hash{$var}{$sample}{clipSeq}{$cPos}; }
                    }
                }
            }
        }
        my $winLen = scalar keys %clippedPos;
        my @sorted = sort { $clippedPos{$b} <=> $clippedPos{$a} } keys %clippedPos;

        if ( $winLen == 0 ) { next; }
        elsif ( $winLen == 1 ) {
            $infile_hash{$var}{leftBkpt}  = $sorted[0];
            $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
        }
        elsif ( $winLen == 2 ) {
            $infile_hash{$var}{leftBkpt}  = $sorted[0];
            $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
            if ( $sorted[1] > $sorted[0] ) {
                $infile_hash{$var}{rightBkpt} = $sorted[1];
            }
            else {
                $infile_hash{$var}{leftBkpt} = $sorted[1];
                $infile_hash{$var}{rightBkpt}--;
            }
        }
        else {
            my @sorted      = sort { $clippedPos{$b} <=> $clippedPos{$a} } keys %clippedPos;
            my $meanClipped = $sumClipped / $winLen;
            my $sumSquares  = 0;
            print "\tmean clipped: $meanClipped\n" if $opts{verbose};

            #for ( my $p = $winStart ; $p <= $winEnd ; $p++ ) {
            #    if ( !defined( $clippedPos{$p} ) ) { $clippedPos{$p} = 0; }
            #    $sumSquares += ( $clippedPos{$p} - $meanClipped )**2;
            #}
            foreach my $cPos ( keys %clippedPos ) {
                $sumSquares += ( $clippedPos{$cPos} - $meanClipped )**2;
            }
            my $sdClipped = sqrt( $sumSquares / ( $winLen - 1 ) );
            print "\tsd clipped: $sdClipped\n"                 if $opts{verbose};
            print "\t1st clipped: $clippedPos{ $sorted[0] }\n" if defined( $sorted[0] ) && $opts{verbose};
            print "\t2nd clipped: $clippedPos{ $sorted[1] }\n" if defined( $sorted[1] ) && $opts{verbose};

            if ( defined( $sorted[0] ) && $clippedPos{ $sorted[0] } > $meanClipped + 1 * $sdClipped ) {
                print "\t*updating to 1st\n" if $opts{verbose};
                $infile_hash{$var}{leftBkpt}  = $sorted[0];
                $infile_hash{$var}{rightBkpt} = $sorted[0] + 1;
            }
            if ( defined( $sorted[1] ) && $clippedPos{ $sorted[1] } > $meanClipped + 1 * $sdClipped ) {
                print "\t*updating to 1st and 2nd\n" if $opts{verbose};
                if ( $sorted[1] > $sorted[0] ) {
                    $infile_hash{$var}{rightBkpt} = $sorted[1];
                }
                else {
                    $infile_hash{$var}{leftBkpt} = $sorted[1];
                    $infile_hash{$var}{rightBkpt}--;
                }
            }
        }
        print "-variant $var completed\n" if $opts{verbose};
    }
    print "Exiting assessBreaks()\n\n" if $opts{verbose};
}

sub scoreData {
    my ( $data_hash, $sample_hash ) = @_;
    print "Entering scoreData()\n" if $opts{verbose};

    foreach my $var ( keys %$data_hash ) {
        print "Variant: $var\n" if $opts{verbose};

        my %priors   = ();
        my %genoFreq = ();
        foreach my $geno ( 0 .. $opts{ploidy} ) {

            #start with uniform priors
            $priors{$geno} = 1 / ( $opts{ploidy} + 1 );
            $genoFreq{$geno}{old} = 0;
            $genoFreq{$geno}{new} = 0;
        }

        my $numIteration = 0;
        my $sumGenoFreq  = 0;
        while ( $numIteration < 10 ) {
            $numIteration++;
            foreach my $sample ( keys %{$sample_hash} ) {

                print "\tSample: $sample\n" if $opts{verbose};
                
                #zero out alternative supporting evidence if below threshold
                if ($$data_hash{$var}{$sample}{numAltRP} + $$data_hash{$var}{$sample}{numAltSR} < $opts{min_evidence}) {
                    $$data_hash{$var}{$sample}{numAltRP} = 0;
                    $$data_hash{$var}{$sample}{numAltSR} = 0;
                    $$data_hash{$var}{$sample}{qualAltRP} = ();
                    $$data_hash{$var}{$sample}{qualAltSR} = ();
                }
                my $numRefRP = $$data_hash{$var}{$sample}{numRefRP};
                my $numAltRP = $$data_hash{$var}{$sample}{numAltRP};
                my $numRefSR = $$data_hash{$var}{$sample}{numRefSR};
                my $numAltSR = $$data_hash{$var}{$sample}{numAltSR};

                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    $$data_hash{$var}{$sample}{pl}{$geno}  = 0;
                    $$data_hash{$var}{$sample}{gl}{$geno}  = 0;
                    $$data_hash{$var}{$sample}{gl0}{$geno} = 0;
                }
                $$data_hash{$var}{$sample}{gq} = 0;
                $$data_hash{$var}{$sample}{gt} = "./.";
                $$data_hash{$var}{$sample}{ft} = "LowQual";

                if ( $numAltRP + $numRefRP + $numAltSR + $numRefSR == 0 ) { next; }
                if ( $$data_hash{$var}{$sample}{avgQ} <= 0 || $$data_hash{$var}{$sample}{avgQ} > 1 ) { $$data_hash{$var}{$sample}{avgQ} = 0.999999; }
                $$data_hash{$var}{$sample}{avgQ} = 0.00001;
                foreach my $g ( 0 .. $opts{ploidy} ) {
                    my $geno = $opts{ploidy} - $g;    #need to reverse as calculation is reference allele based
                    if ( $numAltRP + $numRefRP > 0 &&  1 / $opts{ploidy} ** ($numAltRP + $numRefRP) > 0) {
                        $$data_hash{$var}{$sample}{gl0}{$geno} += calcGl( $opts{ploidy}, $g, $numAltRP + $numRefRP, $numRefRP, $$data_hash{$var}{$sample}{qualRefRP}, $$data_hash{$var}{$sample}{qualAltRP} );
                    }
                    if ( $numAltSR + $numRefSR > 0 && $opts{breakpoint} &&  1 / $opts{ploidy} ** ($numAltSR + $numRefSR) > 0) {
                        $$data_hash{$var}{$sample}{gl0}{$geno} += calcGl( $opts{ploidy}, $g, $numAltSR + $numRefSR, $numRefSR, $$data_hash{$var}{$sample}{qualRefSR}, $$data_hash{$var}{$sample}{qualAltSR} );
                    }
                    print "\tgl0 returned from calcGL() for geno $geno: $$data_hash{$var}{$sample}{gl0}{$geno}\n" if $opts{verbose};
                    if ( $$data_hash{$var}{$sample}{gl0}{$geno} < -255 ) { $$data_hash{$var}{$sample}{gl0}{$geno} = -255; }    #capped
                    $$data_hash{$var}{$sample}{gl}{$geno} = log10( $priors{$geno} ) + $$data_hash{$var}{$sample}{gl0}{$geno};  #adjust with population inferred prior
                }

                my @sortedGeno = sort { $$data_hash{$var}{$sample}{gl}{$b} <=> $$data_hash{$var}{$sample}{gl}{$a} } keys %{ $$data_hash{$var}{$sample}{gl} };

                #calculate PL from GL
                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    print "\tcalculating pl from geno ($$data_hash{$var}{$sample}{gl}{$geno})\n" if $opts{verbose};
                    $$data_hash{$var}{$sample}{pl}{$geno} = int( -10 * $$data_hash{$var}{$sample}{gl}{$geno} );
                    if ( $$data_hash{$var}{$sample}{pl}{$geno} > 255 ) { $$data_hash{$var}{$sample}{pl}{$geno} = 255; }
                    print "\t...$$data_hash{$var}{$sample}{pl}{$geno}\n" if $opts{verbose};
                }

                #normalize PL to most likely genotype
                foreach my $geno ( 0 .. $opts{ploidy} ) {
                    $$data_hash{$var}{$sample}{pl}{$geno} -= $$data_hash{$var}{$sample}{pl}{ $sortedGeno[0] };
                }

                #determine genotype quality
                $$data_hash{$var}{$sample}{gq} = int( 10 * ( $$data_hash{$var}{$sample}{gl}{ $sortedGeno[0] } - $$data_hash{$var}{$sample}{gl}{ $sortedGeno[1] } ) );
                print "\t...$$data_hash{$var}{$sample}{gq}\n" if $opts{verbose};

                my $gt = "0/0";
                if    ( $sortedGeno[0] == 1 ) { $gt = "0/1"; }
                elsif ( $sortedGeno[0] == 2 ) { $gt = "1/1"; }
                $$data_hash{$var}{$sample}{gt} = $gt;
                if   ( $numAltRP + $numRefRP + $numAltSR + $numRefSR < $opts{filter_depth} || $$data_hash{$var}{$sample}{gq} < $opts{filter_qual} ) { $$data_hash{$var}{$sample}{ft} = "LowQual"; }
                else                                                                                                                                { $$data_hash{$var}{$sample}{ft} = "PASS"; }

                $genoFreq{ $sortedGeno[0] }{new}++;
                $sumGenoFreq++;
            }

            my $isConverged = 1;
            if ( $sumGenoFreq == 0 ) { warn "no genotypes passing filters\n"; last; }
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                if ( $genoFreq{$geno}{new} == 0 ) {
                    $priors{$geno} = 1 / $sumGenoFreq;
                }
                else {
                    $priors{$geno} = $genoFreq{$geno}{new} / $sumGenoFreq;
                }
                if ( $genoFreq{$geno}{old} != $genoFreq{$geno}{new} ) { $isConverged = 0; }
                $genoFreq{$geno}{old} = $genoFreq{$geno}{new};
                $genoFreq{$geno}{new} = 0;
            }
            if ($isConverged || !$opts{use_priors}) { last; }
        }
    }
    print "Exiting scoreData()\n\n" if $opts{verbose};
}

sub getDate {
    my ( $second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings ) = localtime();
    my $year = 1900 + $yearOffset;
    $month++;
    my $fmonth = sprintf( "%.2d", $month );
    my $fday   = sprintf( "%.2d", $dayOfMonth );
    return "$year$fmonth$fday";
}

sub report {
    my ( $infile_hash, $data_hash ) = @_;
    print "Entering report()\n" if $opts{verbose};

    #open output file
    if ( defined( $opts{output_filename} ) ) {
        open( foutname1, ">$opts{output_filename}" ) or die("error opening file $opts{output_filename}\n");
    }
    else {
        open( foutname1, ">&", \*STDOUT ) or die;
    }

    my $filedate = getDate();
    print foutname1 <<HEADER;
##fileformat=VCFv4.1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS:MT,Description="Nuclear Mitochondrial Insertion">
##FILTER=<ID=LowQual,Description="No PASS calls in any sample after merging">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number likelihoods">
##FORMAT=<ID=CNL0,Number=.,Type=Float,Description="Copy number likelihoods with no frequency prior">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL0,Number=G,Type=Float,Description="Genotype likelihoods with no frequency prior">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=MSTART,Number=1,Type=Flag,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=MEND,Number=1,Type=Flag,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Flag,Description="Estimated length of mitochondrial insert">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample(s) in which site was originally discovered">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##contig=<ID=GL000207.1,length=4262,assembly=b37>
##contig=<ID=GL000226.1,length=15008,assembly=b37>
##contig=<ID=GL000229.1,length=19913,assembly=b37>
##contig=<ID=GL000231.1,length=27386,assembly=b37>
##contig=<ID=GL000210.1,length=27682,assembly=b37>
##contig=<ID=GL000239.1,length=33824,assembly=b37>
##contig=<ID=GL000235.1,length=34474,assembly=b37>
##contig=<ID=GL000201.1,length=36148,assembly=b37>
##contig=<ID=GL000247.1,length=36422,assembly=b37>
##contig=<ID=GL000245.1,length=36651,assembly=b37>
##contig=<ID=GL000197.1,length=37175,assembly=b37>
##contig=<ID=GL000203.1,length=37498,assembly=b37>
##contig=<ID=GL000246.1,length=38154,assembly=b37>
##contig=<ID=GL000249.1,length=38502,assembly=b37>
##contig=<ID=GL000196.1,length=38914,assembly=b37>
##contig=<ID=GL000248.1,length=39786,assembly=b37>
##contig=<ID=GL000244.1,length=39929,assembly=b37>
##contig=<ID=GL000238.1,length=39939,assembly=b37>
##contig=<ID=GL000202.1,length=40103,assembly=b37>
##contig=<ID=GL000234.1,length=40531,assembly=b37>
##contig=<ID=GL000232.1,length=40652,assembly=b37>
##contig=<ID=GL000206.1,length=41001,assembly=b37>
##contig=<ID=GL000240.1,length=41933,assembly=b37>
##contig=<ID=GL000236.1,length=41934,assembly=b37>
##contig=<ID=GL000241.1,length=42152,assembly=b37>
##contig=<ID=GL000243.1,length=43341,assembly=b37>
##contig=<ID=GL000242.1,length=43523,assembly=b37>
##contig=<ID=GL000230.1,length=43691,assembly=b37>
##contig=<ID=GL000237.1,length=45867,assembly=b37>
##contig=<ID=GL000233.1,length=45941,assembly=b37>
##contig=<ID=GL000204.1,length=81310,assembly=b37>
##contig=<ID=GL000198.1,length=90085,assembly=b37>
##contig=<ID=GL000208.1,length=92689,assembly=b37>
##contig=<ID=GL000191.1,length=106433,assembly=b37>
##contig=<ID=GL000227.1,length=128374,assembly=b37>
##contig=<ID=GL000228.1,length=129120,assembly=b37>
##contig=<ID=GL000214.1,length=137718,assembly=b37>
##contig=<ID=GL000221.1,length=155397,assembly=b37>
##contig=<ID=GL000209.1,length=159169,assembly=b37>
##contig=<ID=GL000218.1,length=161147,assembly=b37>
##contig=<ID=GL000220.1,length=161802,assembly=b37>
##contig=<ID=GL000213.1,length=164239,assembly=b37>
##contig=<ID=GL000211.1,length=166566,assembly=b37>
##contig=<ID=GL000199.1,length=169874,assembly=b37>
##contig=<ID=GL000217.1,length=172149,assembly=b37>
##contig=<ID=GL000216.1,length=172294,assembly=b37>
##contig=<ID=GL000215.1,length=172545,assembly=b37>
##contig=<ID=GL000205.1,length=174588,assembly=b37>
##contig=<ID=GL000219.1,length=179198,assembly=b37>
##contig=<ID=GL000224.1,length=179693,assembly=b37>
##contig=<ID=GL000223.1,length=180455,assembly=b37>
##contig=<ID=GL000195.1,length=182896,assembly=b37>
##contig=<ID=GL000212.1,length=186858,assembly=b37>
##contig=<ID=GL000222.1,length=186861,assembly=b37>
##contig=<ID=GL000200.1,length=187035,assembly=b37>
##contig=<ID=GL000193.1,length=189789,assembly=b37>
##contig=<ID=GL000194.1,length=191469,assembly=b37>
##contig=<ID=GL000225.1,length=211173,assembly=b37>
##contig=<ID=GL000192.1,length=547496,assembly=b37>
##contig=<ID=NC_007605,length=171823,assembly=b37>
##contig=<ID=hs37d5,length=35477943,assembly=b37>
##fileDate=$filedate
##reference=hs37d5
##source=gnomit-$version
HEADER

    my @vars = sort { $$infile_hash{$a}{chr} cmp $$infile_hash{$b}{chr} || $$infile_hash{$a}{start} <=> $$infile_hash{$b}{start} } keys %$infile_hash;
    my @samples = sort { $a cmp $b } keys %{ $data_hash{ $vars[0] } };
    print foutname1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    foreach my $sample (@samples) {
        print foutname1 "\t$sample";
    }
    print foutname1 "\n";
    my $index = 1;
    foreach my $var (@vars) {
        my $chrom  = $$infile_hash{$var}{chr};
        my $id     = $$infile_hash{$var}{id};
        my $alt    = "<INS:MT>";
        my $qual   = $$infile_hash{$var}{qual};
        my $filter = $$infile_hash{$var}{filter};
        my %info   = ();
        my $pos;
        my $end;
        my $ciDelta;

        if ( defined( $$infile_hash{$var}{leftBkpt} ) ) {
            $pos     = $$infile_hash{$var}{leftBkpt} - 1;
            $end     = $$infile_hash{$var}{rightBkpt};
            $ciDelta = $end - $pos + 1;
        }
        else {
            $pos             = $$infile_hash{$var}{start} - 1;
            $end             = $$infile_hash{$var}{end};
            $ciDelta         = $end - $pos + 1;
            $info{IMPRECISE} = undef;
        }
        if ( defined($$infile_hash{$var}{mlen}) ){
            $info{MLEN} = $$infile_hash{$var}{mlen};
            $info{MSTART} = $$infile_hash{$var}{mstart};
            $info{MEND} = $$infile_hash{$var}{mend};
        }
        $info{SAMPLES} = $$infile_hash{$var}{samples};
        $info{CIPOS}  = "0,$ciDelta";
        $info{CIEND}  = "-$ciDelta,0";
        $info{END}    = $end;
        $info{SVTYPE} = "INS";

        my $refline = `$opts{samtools} faidx $opts{reference} $chrom:$pos-$pos`;
        my $ref = ( split( /\n/, $refline ) )[1];
        if ( !defined($ref) ) { $ref = "N"; }
        my $format = "GT:FT:GL0:GQ:PL";

        my $info = "";
        my @sKeys = sort { $a cmp $b } keys %info;
        for ( my $i = 0 ; $i <= $#sKeys ; $i++ ) {
            if ( $i > 0 ) { $info .= ";"; }
            if ( defined( $info{ $sKeys[$i] } ) ) {
                $info .= "$sKeys[$i]=$info{$sKeys[$i]}";
            }
            else {
                $info .= "$sKeys[$i]";
            }
        }
        print foutname1 "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";
        foreach my $sample (@samples) {
            my @gls = ();
            my @pls = ();
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                push @gls, sprintf( "%.2f", $$data_hash{$var}{$sample}{gl0}{$geno} );
                push @pls, $$data_hash{$var}{$sample}{pl}{$geno};
            }
            my $gl = join( ",", @gls );
            my $pl = join( ",", @pls );
            print foutname1 "\t$$data_hash{$var}{$sample}{gt}:$$data_hash{$var}{$sample}{ft}:$gl:$$data_hash{$var}{$sample}{gq}:$pl";
        }
        print foutname1 "\n";
        $index++;
    }
    close(foutname1);
    print "Exiting report()\n\n" if $opts{verbose};
}

sub calcGl {
    my ( $m, $g, $k, $l, $er, $ea ) = @_;
    print "in calcGl():\n"               if $opts{verbose};
    print "\t$m\t$g\t$k\t$l\t$er\t$ea\n" if $opts{verbose};
    if ( 1 / $m**$k <= 0 ) { die "problem in calcGL 1, \t$m\t$g\t$k\t$l\t$er\t$ea\n"; }
    my $gl = log10( 1 / ( $m**$k ) );
    foreach my $e ( @{$er} ) {
        if ( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) <= 0 ) { die "problem in calcGL 2, \t$m\t$g\t$k\t$l\t$e\n"; }
        $gl += log10( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) );
    }
    foreach my $e ( @{$ea} ) {
        if ( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) <= 0 ) { die "problem in calcGL 3, \t$m\t$g\t$k\t$l\t$e\n"; }
        $gl += log10( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) );
    }
    return $gl;
}

sub log10 {
    my $n = shift;
    return log($n) / log(10);
}

sub getInput {
    my ( $infile_hash, $mask_hash, $sample_hash ) = @_;

    print "Entering getInput()\n" if $opts{verbose};

    #input sample information
    open( INFO, $opts{info_filename} ) || die "error opening file $opts{info_filename}, $!\n";
    my @header = ();
    my %fields = ();
    while (<INFO>) {
        chomp;
        if (/^sample/) { @header = split(/\t/); next; }
        my @row = split(/\t/);
        for ( my $i = 0 ; $i <= $#header ; $i++ ) {
            $fields{ $header[$i] } = $row[$i];
        }
        my $sample = $fields{sample};
        $$sample_hash{$sample}{pop}      = $fields{pop};
        $$sample_hash{$sample}{filename} = $fields{filename};
        if ( $fields{median_insert_size} ne "NA" ) {
            $$sample_hash{$sample}{winlen}                    = $fields{median_insert_size} + 3 * $fields{median_absolute_deviation};
            $$sample_hash{$sample}{median_insert_size}        = $fields{median_insert_size};
            $$sample_hash{$sample}{median_absolute_deviation} = $fields{median_absolute_deviation};
        }
        else {
            $$sample_hash{$sample}{winlen}                    = 500;
            $$sample_hash{$sample}{median_insert_size}        = 400;
            $$sample_hash{$sample}{median_absolute_deviation} = 40;
        }
        $$sample_hash{$sample}{meancoverage} = $fields{mean_coverage};
        if ( defined( $fields{read_groups} ) ) {
            my @rgs = split( /,/, $fields{read_groups} );
            %{ $$sample_hash{$sample}{read_groups} } = map { $_, 1 } @rgs;
        }
    }
    close INFO;

    #input mask coordinates
    open( MASK, $opts{mask_filename} ) || die "Could not open $opts{mask_filename} for input, $!\n";
    while (<MASK>) {
        chomp;
        my ( $chr, $start, $end, $id ) = split(/\t/);
        $chr =~ s/chr//g;
        $$mask_hash{$chr}{$start} = $end;
    }
    close MASK;

    #input variant coordinates
    if ( $opts{input_filename} =~ /\.gz$/ ) {
        open( VARS, "zcat $opts{input_filename} |" ) || die "Could not open $opts{input_filename} for input, $!\n";
    }
    else {
        open( VARS, $opts{input_filename} ) || die "Could not open $opts{input_filename} for input, $!\n";
    }
    my $varnum = 1;

    #BED FORMAT
    #while (<VARS>) {
    #    chomp;
    #    my ( $chr, $start, $end ) = split(/\t/);
    #    $chr =~ s/chr//g;
    #    $$infile_hash{$varnum}{chr}   = $chr;
    #    $$infile_hash{$varnum}{start} = $start;
    #    $$infile_hash{$varnum}{end}   = $end;
    #    $varnum++;
    #}

    #VCF FORMAT
    while (<VARS>) {
        next if /^#/;
        chomp;
        my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info ) = split(/\t/);
        $chr =~ s/chr//g;
        if ( defined( $opts{chr} ) ) {
            if ( $chr ne $opts{chr} ) { next; }
        }
        $infile_hash{$varnum}{chr}   = $chr;
        $infile_hash{$varnum}{id}    = $id;
        $infile_hash{$varnum}{start} = $pos + 1;
        my ($end) = $info =~ /END=(\d+)/;
        my ($mstart) = $info =~ /MSTART=(\d+)/;
        my ($mend) = $info =~ /MEND=(\d+)/;
        my ($mlen) = $info =~ /MLEN=(\d+)/;
        my ($samples) = $info =~ /SAMPLES=(\w+?);/;
        if (defined($mlen)) {
            $infile_hash{$varnum}{mlen} = $mlen;
            $infile_hash{$varnum}{mstart} = $mstart;
            $infile_hash{$varnum}{mend} = $mend;
        }
        $infile_hash{$varnum}{samples} = $samples;
        $infile_hash{$varnum}{end}    = $end;
        $infile_hash{$varnum}{filter} = $filter;
        $infile_hash{$varnum}{qual}   = $qual;
        $varnum++;
    }
    close VARS;
    print "Exiting getInput()\n\n" if $opts{verbose};
}

sub getMateInfo {
    my ( $qname, $rnext, $pnext, $readgroup_hash, $filename ) = @_;

    my $command = "";
    if ( $opts{by_chr_dir} ) {
        if ( $opts{ucsc} ) {
            $command = "$opts{samtools} view $filename" . "chr$rnext.*bam chr$rnext:$pnext-$pnext |";
        }
        else {
            $command = "$opts{samtools} view $filename$rnext.*bam $rnext:$pnext-$pnext |";
        }
    }
    else {
        if ( $opts{ucsc} ) {
            $command = "$opts{samtools} view $filename chr$rnext:$pnext-$pnext |";
        }
        else {
            $command = "$opts{samtools} view $filename $rnext:$pnext-$pnext |";
        }
    }
    open( MCI, "$command" ) || die "Could not open $filename, $!\n";

    my $cFlag    = 0;
    my $cPos     = -1;
    my $clipside = "n";
    my $clipsize = -1;
    my $clipseq = "";
    my $seqLen   = 0;
    my $seq      = "";
    my $matchLen = 0;

    while (<MCI>) {
        chomp;
        my ( $m_qname, $m_flag, $m_rname, $m_pos, $m_mapq, $m_cigar, $m_rnext, $m_pnext, $m_tlen, $m_seq, $m_qual, $opt ) = split(/\t/);
        if ( $m_qname ne $qname ) { next; }
        my ($read_group) = $_ =~ /RG:Z:(\S+)/;

        if ( $opts{read_groups} && !defined($read_group) ) { next; }
        elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

        ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $m_pos, $m_cigar, $m_qual, $m_seq );
        $seqLen   = 0;
        $matchLen = 0;
        $seq      = $m_seq;
        while ( $m_cigar =~ /(\d+)M/g ) {
            $seqLen   += $1;
            $matchLen += $1;
        }
        while ( $m_cigar =~ /(\d+)N/g ) {
            $seqLen += $1;
        }
        while ( $m_cigar =~ /(\d+)D/g ) {
            $seqLen += $1;
        }
    }
    close MCI;

    return ( $cFlag, $cPos, $clipside, $clipsize, $seqLen, $matchLen, $seq );
}

sub refineData {
    my ( $infile_hash, $mask_hash, $sample_hash, $data_hash ) = @_;
    print "Entering refineData()\n" if $opts{verbose};

    my $numSamp = 0;
    foreach my $var ( keys %{$infile_hash} ) {
        if ( !defined( $$infile_hash{$var}{leftBkpt} ) ) { next; }
        foreach my $sample ( keys %{$sample_hash} ) {
            $numSamp++;

            #last if $numSamp > 50;

            #currently basing coordinates off of ONLY start position, may need to revisit
            my $chr     = $$infile_hash{$var}{chr};
            my $l_start = $$infile_hash{$var}{leftBkpt} - $$sample_hash{$sample}{winlen};
            my $l_end   = $$infile_hash{$var}{leftBkpt};
            my $r_start = $$infile_hash{$var}{rightBkpt};
            my $r_end   = $$infile_hash{$var}{rightBkpt} + $$sample_hash{$sample}{winlen};
            my $command = "";
            my $chrM    = "";

            print "REFINED: $var\t$chr\t$l_start\t$l_end\t$r_start\t$r_end\n" if $opts{verbose};
            if ( $opts{by_chr_dir} ) {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename}chr$chr.*bam chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename}$chr.*bam $chr:$l_start-$r_end |";
                    $chrM    = "MT";
                }
            }
            else {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename} chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename} $chr:$l_start-$r_end |";
                    $chrM    = "MT";
                }
            }
            $$data_hash{$var}{$sample}{numRefRP}  = 0;
            $$data_hash{$var}{$sample}{numAltRP}  = 0;
            $$data_hash{$var}{$sample}{numRefSR}  = 0;
            $$data_hash{$var}{$sample}{numAltSR}  = 0;
            $$data_hash{$var}{$sample}{qualRefRP} = ();
            $$data_hash{$var}{$sample}{qualAltRP} = ();
            $$data_hash{$var}{$sample}{qualRefSR} = ();
            $$data_hash{$var}{$sample}{qualAltSR} = ();
            $$data_hash{$var}{$sample}{avgQ}      = 0;

            my %found = ();

            #print "command: $command\n" if $opts{verbose};
            open( SAM, $command ) || die "error in opening file, $!\n";
            while (<SAM>) {
                chomp;
                my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);
                $rname =~ s/chr//g;
                $rnext =~ s/chr//g;

                my $mapE = 10**( -1 * $mapq / 10 );
                my ($read_group) = $_ =~ /RG:Z:(\S+)/;

                if ( $opts{read_groups} && !defined($read_group) ) { next; }
                elsif ( $opts{read_groups} && !defined( $$sample_hash{$sample}{read_groups}{$read_group} ) ) { next; }

                if ( $mapq < $opts{min_map_qual} ) { next; }

                my $dir = 0;    #F
                if ( $flag & 16 ) { $dir = 1; }    #R
                my $dnext = 0;                     #mate F
                if ( $flag & 32 ) { $dnext = 1; }  #mate R

                my $seqLen   = 0;
                my $matchLen = 0;
                while ( $cigar =~ /(\d+)M/g ) {
                    $matchLen += $1;
                    $seqLen   += $1;
                }
                while ( $cigar =~ /(\d+)N/g ) {
                    $seqLen += $1;
                }
                while ( $cigar =~ /(\d+)D/g ) {
                    $seqLen += $1;
                }

                my ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $pos, $cigar, $qual, $seq );
                my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize, $m_seqLen, $m_matchLen, $m_seq ) = getMateInfo( $qname, $rname, $pnext, $$sample_hash{$sample}{read_groups}, $$sample_hash{$sample}{filename} );

                if ( $cFlag == 1 && ( abs( $cPos - $$infile_hash{$var}{leftBkpt} ) <= $opts{min_clipped_seq} || abs( $cPos - $$infile_hash{$var}{rightBkpt} ) <= $opts{min_clipped_seq} ) ) {

                    $$data_hash{$var}{$sample}{numAltSR}++;
                    push @{ $$data_hash{$var}{$sample}{qualAltSR} }, $mapE;
                }
                elsif ( $matchLen >= 0.95 * length($seq) && ( ( $pos < $$infile_hash{$var}{leftBkpt} && $pos + $seqLen - 1 > $$infile_hash{$var}{leftBkpt} && abs( $$infile_hash{$var}{leftBkpt} - $pos ) > $opts{min_clipped_seq} && abs( $$infile_hash{$var}{leftBkpt} - ( $pos + $seqLen - 1 ) ) > $opts{min_clipped_seq} ) || ( $pos < $$infile_hash{$var}{rightBkpt} && $pos + $seqLen - 1 > $$infile_hash{$var}{rightBkpt} && abs( $$infile_hash{$var}{rightBkpt} - $pos ) > $opts{min_clipped_seq} && abs( $$infile_hash{$var}{rightBkpt} - ( $pos + $seqLen - 1 ) ) > $opts{min_clipped_seq} ) ) ) {

                    #non-clipped reads with slight overlap with breakpoint may have mismatches instead of clipping, so skip
                    #also require at least 95% of read be 'matched' to reference, as the presence of insertion can cause some wacky mappings
                    $$data_hash{$var}{$sample}{numRefSR}++;
                    push @{ $$data_hash{$var}{$sample}{qualRefSR} }, $mapE;
                }

                if ( $rnext eq $chrM || checkMaskOverlap( $chr, $pos, $mask_hash ) ) {

                    #mate mapped to mt sequence
                    if ( ( $dir == 0 && $pos >= $l_start && $pos <= $l_end ) || ( $dir == 1 && $pos >= $r_start && $pos <= $r_end ) ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                }
                elsif ( $dir == 0 && $dnext == 1 && $rnext eq "=" && !$found{$qname} ) {
                    $found{$qname} = 1;
                    if ( abs($tlen) > $$sample_hash{$sample}{median_insert_size} + 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of upper bounds
                    if ( abs($tlen) < $$sample_hash{$sample}{median_insert_size} - 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of lower bounds

                    #allow some flexibility for clipped positions around breakpoint due to alignment artifacts (later may want to inplement a realignment step)
                    if    ( $cPos > -1   && abs( $cPos - $l_end ) <= $opts{min_clipped_seq}     && $cFlag == 1 )   { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $m_cPos > -1 && abs( $m_cPos - $l_end ) <= $opts{min_clipped_seq}   && $m_cFlag == 1 ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $cPos > -1   && abs( $cPos - $r_start ) <= $opts{min_clipped_seq}   && $cFlag == 1 )   { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( $m_cPos > -1 && abs( $m_cPos - $r_start ) <= $opts{min_clipped_seq} && $m_cFlag == 1 ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                    elsif ( ( $pos <= $l_end && $pnext + $m_seqLen >= $l_end ) || ( $pos <= $r_start && $pnext + $m_seqLen >= $r_start ) ) {

                        #non-clipped reads with slight overlap with breakpoint may have mismatches instead of clipping
                        if (   abs( $pos - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pos + $seqLen - 1 - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pnext - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pnext + $m_seqLen - 1 - $l_end ) > $opts{min_clipped_seq}
                            && abs( $pos - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pos + $seqLen - 1 - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pnext - $r_start ) > $opts{min_clipped_seq}
                            && abs( $pnext + $m_seqLen - 1 - $r_start ) > $opts{min_clipped_seq}
                            && $matchLen >= 0.95 * length($seq)
                            && $m_matchLen >= 0.95 * length($m_seq) )
                        {
                            $data_hash{$var}{$sample}{numRefRP}++;
                            push @{ $$data_hash{$var}{$sample}{qualRefRP} }, $mapE;
                        }
                    }
                }
                $$data_hash{$var}{$sample}{avgQ} += $mapE;
            }
            close SAM;

            if ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} > 0 ) {
                $$data_hash{$var}{$sample}{avgQ} /= ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} );
            }
            else {
                $data_hash{$var}{$sample}{avgQ} = 0.00001;
            }
            print "$sample\t$var\t$$data_hash{$var}{$sample}{numRefRP}\t$$data_hash{$var}{$sample}{numRefSR}\t$$data_hash{$var}{$sample}{numAltRP}\t$$data_hash{$var}{$sample}{numAltSR}\t$$data_hash{$var}{$sample}{avgQ}\n" if $opts{verbose};
        }
    }
    print "Exiting refineData()\n\n" if $opts{verbose};
}

sub getData {
    my ( $infile_hash, $mask_hash, $sample_hash, $data_hash ) = @_;
    print "Entering getData()\n" if $opts{verbose};

    my $numSamp = 0;
    foreach my $var ( keys %{$infile_hash} ) {
        foreach my $sample ( keys %{$sample_hash} ) {
            $numSamp++;
            my $chr     = $$infile_hash{$var}{chr};
            my $l_start = $$infile_hash{$var}{start} - $$sample_hash{$sample}{winlen};
            my $l_end   = $$infile_hash{$var}{start};
            my $r_start = $$infile_hash{$var}{end};
            my $r_end   = $$infile_hash{$var}{end} + $$sample_hash{$sample}{winlen};
            my $command = "";
            my $chrM    = "";

            #print "$var\t$chr\t$l_start\t$l_end\t$r_start\t$r_end\n" if $opts{verbose};
            if ( $opts{by_chr_dir} ) {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename}chr$chr.*bam chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename}$chr.*bam $chr:$l_start-$r_end |";
                    $chrM    = "MT";
                }
            }
            else {
                if ( $opts{ucsc} ) {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename} chr$chr:$l_start-$r_end |";
                    $chrM    = "M";
                }
                else {
                    $command = "$opts{samtools} view $$sample_hash{$sample}{filename} $chr:$l_start-$r_end |";
                    $chrM    = "MT";
                }
            }
            $$data_hash{$var}{$sample}{numRefRP}  = 0;
            $$data_hash{$var}{$sample}{numAltRP}  = 0;
            $$data_hash{$var}{$sample}{numRefSR}  = 0;
            $$data_hash{$var}{$sample}{numAltSR}  = 0;
            $$data_hash{$var}{$sample}{qualRefRP} = ();
            $$data_hash{$var}{$sample}{qualAltRP} = ();
            $$data_hash{$var}{$sample}{qualRefSR} = ();
            $$data_hash{$var}{$sample}{qualAltSR} = ();
            $$data_hash{$var}{$sample}{avgQ}      = 0.000001;

            #print "command: $command\n" if $opts{verbose};
            open( SAM, $command ) || die "error in opening file, $!\n";
            while (<SAM>) {
                chomp;
                my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);
                $rname =~ s/chr//g;
                $rnext =~ s/chr//g;

                my $mapE = 10**( -1 * $mapq / 10 );
                my ($read_group) = $_ =~ /RG:Z:(\S+)/;

                if ( $opts{read_groups} && !defined($read_group) ) { next; }
                elsif ( $opts{read_groups} && !defined( $$sample_hash{$sample}{read_groups}{$read_group} ) ) { next; }

                if ( $mapq < $opts{min_map_qual} ) { next; }

                my $dir = 0;    #F
                if ( $flag & 16 ) { $dir = 1; }    #R
                my $dnext = 0;                     #mate F
                if ( $flag & 32 ) { $dnext = 1; }  #mate R

                my ( $cFlag, $cPos, $clipside, $clipsize, $clipseq ) = getSoftClipInfo( $pos, $cigar, $qual, $seq );
                if ( $cPos > -1 && $cFlag == 1 && $cPos >= $$infile_hash{$var}{start} - $opts{clipped_flank} && $cPos <= $$infile_hash{$var}{end} + $opts{clipped_flank} ) {

                    #potential breakpoints must be at or between previous breakpoint bounds
                    $$data_hash{$var}{$sample}{clipPos}{$cPos}{$cFlag}++;
                    if ( !defined( $$data_hash{$var}{$sample}{clipSeq}{$cPos} ) || length($clipseq) > length($$data_hash{$var}{$sample}{clipSeq}{$cPos}) ) {
                        $$data_hash{$var}{$sample}{clipSeq}{$cPos} = $clipseq;
                    }
                }
                if ( $rnext eq $chrM || checkMaskOverlap( $chr, $pos, $mask_hash ) ) {

                    #mate mapped to mt sequence
                    if ( ( $dir == 0 && $pos >= $l_start && $pos <= $l_end ) || ( $dir == 1 && $pos >= $r_start && $pos <= $r_end ) ) { $$data_hash{$var}{$sample}{numAltRP}++; push @{ $$data_hash{$var}{$sample}{qualAltRP} }, $mapE; }
                }
                elsif ( $dir == 0 && $dnext == 1 && $rnext eq "=" ) {
                    if ( abs($tlen) > $$sample_hash{$sample}{median_insert_size} + 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of upper bounds
                    if ( abs($tlen) < $$sample_hash{$sample}{median_insert_size} - 3 * $$sample_hash{$sample}{median_absolute_deviation} ) { next; }    #insert length of potential reference supporting allele out of lower bounds
                                                                                                                                                        #"normal" reads which clip at the breakpoint are -not- reference supporting!
                                                                                                                                                        # --------->..............<-------L------ or <------L///////
                    if ( $pos >= $l_start && $pos < $l_end && $pnext < $l_end && $pnext + length($seq) > $l_end ) {

                        #my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize ) = getMateInfo( $qname, $rname, $pnext, $readgroup_hash, $$sample_hash{$sample}{filename});
                        #if   ( $m_cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                  { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # --------->..............<-------R------ or <//////R-------
                    elsif ( $pos >= $l_start && $pos < $l_end && $pnext <= $r_end && $pnext + length($seq) > $r_end ) {

                        #my ( $m_cFlag, $m_cPos, $m_clipside, $m_clipsize ) = getMateInfo( $qname, $rname, $pnext, $readgroup_hash, $$sample_hash{$sample}{filename} );
                        #if   ( $m_cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                  { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # --------L//////> or -------L----->..............<----------
                    elsif ( $pos >= $l_start && $pos < $l_end && $pnext > $l_end && $pos + length($seq) > $l_end ) {

                        #if   ( $cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # \\\\\\\\R------> or -------R----->..............<----------
                    elsif ( $pos <= $r_start && $pos + length($seq) > $r_start && $pnext > $r_start ) {

                        #if   ( $cPos > -1 ) { $$data_hash{$var}{$sample}{numAltRP}++; }
                        #else                { $$data_hash{$var}{$sample}{numRefRP}++; }
                    }

                    # ------------->......L......<------------- or ------------->......R......<-------------
                    elsif ( ( $pos >= $l_start && $pos < $l_end && $pnext > $l_end && $pnext <= $r_end ) || ( $pos >= $l_start && $pos < $r_start && $pnext > $r_start && $pnext <= $r_end ) ) {
                        $data_hash{$var}{$sample}{numRefRP}++;
                        push @{ $$data_hash{$var}{$sample}{qualRefRP} }, $mapE;
                    }
                }
                $$data_hash{$var}{$sample}{avgQ} += $mapE;
            }
            close SAM;

            if ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} > 0 ) {
                $$data_hash{$var}{$sample}{avgQ} /= ( $$data_hash{$var}{$sample}{numRefRP} + $$data_hash{$var}{$sample}{numAltRP} );
            }
            else {
                $data_hash{$var}{$sample}{avgQ} = 0.99;
            }
            print "\t$sample\t$var\t$$data_hash{$var}{$sample}{numRefRP}\t$$data_hash{$var}{$sample}{numAltRP}\t$$data_hash{$var}{$sample}{avgQ}\n" if $opts{verbose};
        }
    }
    print "Exiting getData()\n\n" if $opts{verbose};
}

sub checkMaskOverlap {
    my ( $chr, $pos, $mask_hash ) = @_;
    my $isMaskOverlap = 0;
    foreach my $maskStart ( keys %{ $$mask_hash{$chr} } ) {
        my $maskEnd = $$mask_hash{$chr}{$maskStart};
        if ( $pos >= $maskStart && $pos <= $maskEnd ) {
            $isMaskOverlap = 1;
            last;
        }
    }
    return $isMaskOverlap;
}

sub getSoftClipInfo {
    my ( $pos, $cigar, $qual, $seq ) = @_;
    my $clipside = "";
    my $clipsize = 0;
    my $clipseq  = "";
    my $cPos     = -1;
    my $avgQual  = -1;
    my $cFlag    = 1;

    if ( $cigar =~ /^(\d+)S.*M.*?(\d+)S$/ ) {
        if ( $1 > $2 ) {
            $cPos     = $pos;
            $clipside = "l";
            $clipsize = $1;

            #consider breakpoint after leftmost soft clipped fragment
        }
        else {
            $cPos     = $pos - 1;
            $clipside = "r";
            $clipsize = $2;
            while ( $cigar =~ /(\d+)M/g ) {    #have to take into account that a CIGAR may contain multiple M's
                $cPos += $1;
            }
            while ( $cigar =~ /(\d+)I/g ) {
                $cPos -= $1;
            }
            while ( $cigar =~ /(\d+)D/g ) {
                $cPos += $1;
            }
        }
    }

    #upstream soft clip only
    elsif ( $cigar =~ /^(\d+)S.*M/ ) {
        $cPos     = $pos;
        $clipside = "l";
        $clipsize = $1;
    }

    #downstream soft clip only
    elsif ( $cigar =~ /M.*?(\d+)S/ ) {
        $cPos     = $pos - 1;
        $clipside = "r";
        $clipsize = $1;
        while ( $cigar =~ /(\d+)M/g ) {    #have to take into account that a CIGAR may contain multiple M's
            $cPos += $1;
        }
        while ( $cigar =~ /(\d+)I/g ) {
            $cPos -= $1;
        }
        while ( $cigar =~ /(\d+)D/g ) {
            $cPos += $1;
        }
    }

    #Check quality of clipped sequence and alignment to reference
    if ( $cPos > -1 ) {
        my $clippedQuals = "";

        if ( $clipside eq "r" ) {
            $clipseq   = substr( $seq,  length($qual) - $clipsize - 1, $clipsize );
            $clippedQuals = substr( $qual, length($qual) - $clipsize - 1, $clipsize );
        }
        else {
            $clipseq   = substr( $seq,  0, $clipsize );
            $clippedQuals = substr( $qual, 0, $clipsize );
        }

        my $avgQualSum = 0;
        my $avgQualNum = 0;
        foreach my $qual ( split( //, $clippedQuals ) ) {
            $avgQualSum += ord($qual) - 33;
            $avgQualNum++;
        }
        $avgQual = $avgQualSum / $avgQualNum;
    }
    if ( $avgQual < 10 ) { $cFlag = 0; }
    if ( $clipsize < $opts{min_clipped_seq} ) { $cFlag = 0; }
    return ( $cFlag, $cPos, $clipside, $clipsize, $clipseq );
}

sub usage {
    my $version = shift;
    printf("\n");
    printf( "%-9s %s\n", "Program:", "gnomit.pl" );
    printf( "%-9s %s\n", "Version:", "$version" );
    printf("\n");
    printf( "%-9s %s\n", "Usage:", "gnomit.pl [options]" );
    printf("\n");
    printf( "%-9s %-35s %s\n", "Options:", "--input_filename=[filename]",     "Input alignment file in BAM format" );
    printf( "%-9s %-35s %s\n", "",         "--info_filename=[filename]",      "Input file wth per-sample information (required)" );
    printf( "%-9s %-35s %s\n", "",         "--output_filename=[filename]",    "Output file (default stdout)" );
    printf( "%-9s %-35s %s\n", "",         "--mask_filename=[filename]",      "Mask file for reference numts in BED format (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--reference=[filename]",          "Reference file");
    printf( "%-9s %-35s %s\n", "",         "--include_mask",                  "Include aberrant reads mapped to mask regions in clustering" );
    printf( "%-9s %-35s %s\n", "",         "--breakpoint",                    "Include soft clipped reads in likelihood calculation" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_include=[integer]", "Maximum distance to be included in cluster (default 600)" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_link=[integer]",    "Maximum distance to link clusters (default 800)" );
    printf( "%-9s %-35s %s\n", "",         "--min_reads_cluster=[integer]",   "Minimum number of reads to link a cluster (default 1)" );
    printf( "%-9s %-35s %s\n", "",         "--min_evidence=[integer]",        "Minimum evidence to consider an insertion event for genotyping (default 3)" );
    printf( "%-9s %-35s %s\n", "",         "--min_map_qual=[integer]",        "Minimum mapping quality for read consideration (default 10)" );
    printf( "%-9s %-35s %s\n", "",         "--max_read_cov=[integer]",        "Maximum read coverage allowed for breakpoint searching (default 200)" );
    printf( "%-9s %-35s %s\n", "",         "--min_clipped_seq=[integer]",     "Minimum clipped sequence required to consider as putative breakpoint (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--max_num_clipped=[integer]",     "Maximum number of clipped sequences observed before removing from evidence consideration (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--read_groups",                   "Stratify analysis to specified read group(s) indicated in info_filename (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--use_priors",                    "Estimate priors using EM framework" );
    printf( "%-9s %-35s %s\n", "",         "--by_chr_dir",                    "If set, expects to find chr specific BAM files info_filename indicated directory" );
    printf( "%-9s %-35s %s\n", "",         "--prefix=[string]",               "Prepend label in report output" );
    printf( "%-9s %-35s %s\n", "",         "--ucsc",                          "Use UCSC genome formatting (e.g. chrM)" );
    printf("\n");
}

sub checkOptions {
    my $optResult = shift;
    my $opts      = shift;
    my $version   = shift;

    if ( !$optResult || $$opts{help} ) {
        usage($version);
        exit;
    }

    if ( !defined( $$opts{input_filename} ) && !defined( $$opts{by_chr_dir} ) ) {
        print "\n***ERROR***\t--input_filename or --by_chr_dir is required\n";
        usage($version);
        exit;
    }
    elsif ( !defined( $$opts{by_chr_dir} ) && !-e $$opts{input_filename} ) {
        print "\n***ERROR***\t--input_filename does not exist\n";
        usage($version);
        exit;
    }
    elsif ( defined( $$opts{by_chr_dir} ) && !-d $$opts{by_chr_dir} ) {
        print "\n***ERROR***\t--by_chr_dir does not exist\n";
        usage($version);
        exit;
    }
    if ( defined( $$opts{mask_filename} ) && !-e ( $$opts{mask_filename} ) ) {
        print "\n***ERROR***\t--mask_filename does not exist\n";
        usage($version);
        exit;
    }
    if ( !$$opts{include_mask} && !defined( $$opts{mask_filename} ) ) {
        print "\n***ERROR***\t--mask_filename is neccessary with --include_mask option\n";
        usage($version);
        exit;
    }
}

sub log2 {
    my $n = shift;
    return log($n) / log(2);
}
