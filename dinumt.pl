#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $version = "0.0.22";

#version update
# 0.0.22
#   -changed name to "dinumt" (dynumite!)
#   -added option --mt_names for discrete MT identifiers
#
# 0.0.21
#   -added option for ensemble genomes (chrMT)
#   -fix usage of masking when --include-mask isn't present
#   -include reference paremeter in vcf output
#
# 0.0.20
#   -added option to output GL information
#
# 0.0.19
#   -added option to output supporting reads to auxilliary file
#
# 0.0.18
#   -added mito position estimation
#
# 0.0.17
#   -implemented likelihood scoring for events
#   -added quality, evidence and depth filters
#   -implemented VCF format reporting
#
# 0.0.15
#   -bug fixes
#
# 0.0.14
#   -fixed bug with parsing read group information
#   -fixed bug where read group information was not being utilized in findBreakpoint()
#
# 0.0.13
#   -added option for minimum clip size to consider
#   -added option for maximum limit of putative breakpoints
#   -added estimated numt size and mitochondria coordinates to output report
#
# 0.0.12
#   -added option to additionally attempt to cluster reads mapping in known numt regions
#    whose mates map elsewhere
#   -added prefix option for report
#
# 0.0.11
#   -added restriction that mates of linked clusters must be consistent with direct
#    or inverted sequence (incl changing dir to dnext in input_hash)
#
# 0.0.10
#   -added option for maximun read coverage when considering breakpoints
#   -moved mate quality filtering to getInput()
#
# 0.0.9
#   -changed default cluster reads to 1
#   -added option to consider total evidence from read pairs and breakpoints
#   -moved mask file comparison to getInput()
#   -added option for minimum mapping quality
#   -added filter in seqcluster() to remove low quality reads/clusters
#
# 0.0.8
#   -added option to restrict analysis to one or more read groups
#   -added option to use UCSC naming conventions (e.g. chrM instead of MT)
#   -added option to use genomes segregated by chromosome

my %opts = ();

$opts{len_cluster_include} = 600;
$opts{len_cluster_link}    = 800;
$opts{filter_quality}      = 50;
$opts{filter_evidence}     = 4;
$opts{filter_depth}        = 5;
$opts{min_reads_cluster}   = 1;
$opts{min_clipped_seq}     = 5;
$opts{max_num_clipped}     = 5;
$opts{include_mask}        = 0;
$opts{min_evidence}        = 4;
$opts{min_map_qual}        = 10;
$opts{max_read_cov}        = 200;
$opts{mask_filename}       = "numtS.bed";
$opts{reference}           = "genome.fa";
$opts{samtools}            = "samtools";
$opts{prefix}              = "numt";
$opts{len_mt}              = 16596;                                                               #eventually should be read in by BAM header
$opts{ploidy}              = 2;
$opts{output_support}      = 0;
$opts{output_gl}           = 0;
my $optResult = GetOptions(
    "input_filename=s"      => \$opts{input_filename},
    "output_filename=s"     => \$opts{output_filename},
    "mask_filename=s"       => \$opts{mask_filename},
    "support_filename=s"    => \$opts{support_filename},
    "include_mask"          => \$opts{include_mask},
    "output_support"        => \$opts{output_support},
    "len_cluster_include=i" => \$opts{len_cluster_include},
    "len_cluster_link=i"    => \$opts{len_cluster_link},
    "filter_quality=i"      => \$opts{filter_quality},
    "filter_evidence=i"     => \$opts{filter_evidence},
    "filter_depth=i"        => \$opts{filter_evidence},
    "min_reads_cluster=i"   => \$opts{min_reads_cluster},
    "min_evidence=i"        => \$opts{min_evidence},
    "min_clipped_seq=i"     => \$opts{min_clipped_seq},
    "max_num_clipped=i"     => \$opts{max_num_clipped},
    "min_map_qual=i"        => \$opts{min_map_qual},
    "max_read_cov=i"        => \$opts{max_read_cov},
    "mean_read_cov=f"       => \$opts{mean_read_cov},
    "insert_size=s"         => \$opts{insert_size},
    "read_groups=s"         => \$opts{read_groups},
    "mt_names=s"            => \$opts{mt_names},
    "by_chr_dir=s"          => \$opts{by_chr_dir},
    "reference=s"           => \$opts{reference},
    "prefix=s"              => \$opts{prefix},
    "output_gl"             => \$opts{output_gl},
    "ucsc"                  => \$opts{ucsc},
    "ensembl"               => \$opts{ensembl},
    "help"                  => \$opts{help},
    "verbose"               => \$opts{verbose}
);

checkOptions( $optResult, \%opts, $version );

my $seq_num  = 0;
my %seq_hash = ();

my %sorted_hash    = ();
my %readgroup_hash = ();

my $i            = 1;
my %infile_hash  = ();
my %group_hash   = ();
my %outfile_hash = ();
my %mask_hash    = ();
my %mt_hash      = ();

if ( defined( $opts{read_groups} ) ) {
    my @rgs = split( /,/, $opts{read_groups} );
    %readgroup_hash = map { $_, 1 } @rgs;
}

if (defined( $opts{mt_names} ) ) {
    my @mts = split (/,/, $opts{mt_names} );
    %mt_hash = map { $_, 1 } @mts;
}

getInput( \%infile_hash, \%readgroup_hash, \%mask_hash, \%mt_hash );
seqCluster( \%infile_hash );
linkCluster( \%infile_hash );
mapCluster( \%infile_hash, \%outfile_hash, \%readgroup_hash );
findBreakpoint( \%outfile_hash, \%readgroup_hash, \%mask_hash );
scoreData( \%outfile_hash );
report( \%outfile_hash );

################################################################################################################
sub scoreData {
    my ($outfile_hash) = @_;
    print "entering scoreData()\n" if $opts{verbose};
    foreach my $group ( keys %$outfile_hash ) {
        print "Group: $group\n" if $opts{verbose};
        my $sumGP  = 0;
        my $numRef = $$outfile_hash{$group}{numRefRP};
        my $numAlt = $$outfile_hash{$group}{numAltRP};
        if ( $$outfile_hash{$group}{numAltSR} > 0 ) {
            $numRef += $$outfile_hash{$group}{numRefSR};
            $numAlt += $$outfile_hash{$group}{numAltSR};
        }
        print "\t$numRef\t$numAlt\n" if $opts{verbose};

        foreach my $g ( 0 .. $opts{ploidy} ) {
            my $geno = $opts{ploidy} - $g;    #need to reverse as calculation is reference allele based
            if ( $numAlt + $numRef > 0 && 1 / $opts{ploidy}**( $numAlt + $numRef ) > 0 ) {
                $$outfile_hash{$group}{gl}{$geno} = calcGl( $opts{ploidy}, $g, $numAlt + $numRef, $numRef, $$outfile_hash{$group}{avgQ} );
                $$outfile_hash{$group}{gp}{$geno} = 10**$$outfile_hash{$group}{gl}{$geno};
                $sumGP += $$outfile_hash{$group}{gp}{$geno};
                print "\t$geno\t$$outfile_hash{$group}{gl}{$geno}\t$$outfile_hash{$group}{gp}{$geno}\n" if $opts{verbose};
            }
        }
        print "\tsumGP: $sumGP\n" if $opts{verbose};
        if ( $sumGP == 0 ) {
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                $$outfile_hash{$group}{pl}{$geno} = 0;
                $$outfile_hash{$group}{gl}{$geno} = 0;
            }
            $$outfile_hash{$group}{gq} = 0;
            $$outfile_hash{$group}{gt} = "./.";
            $$outfile_hash{$group}{ft} = "NC";
        }
        else {
            my $maxGP   = 0;
            my $maxGeno = 0;
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                if ( $$outfile_hash{$group}{gp}{$geno} == 0 ) { $$outfile_hash{$group}{gp}{$geno} = 1e-200; }
                $$outfile_hash{$group}{gp}{$geno} /= $sumGP;
                $$outfile_hash{$group}{pl}{$geno} = int( -10 * log10( $$outfile_hash{$group}{gp}{$geno} ) );
                if ( $$outfile_hash{$group}{gp}{$geno} > $maxGP ) { $maxGP = $$outfile_hash{$group}{gp}{$geno}; $maxGeno = $geno; }
            }

            $maxGP = 1 - $$outfile_hash{$group}{gp}{0};    #calculate P(not 0/0 | data)
            if ( 1 - $maxGP == 0 ) {
                $$outfile_hash{$group}{gq} = 199;
            }
            else {
                $$outfile_hash{$group}{gq} = int( -10 * log10( 1 - $maxGP ) );
            }
            my $gt = "0/0";
            if    ( $maxGeno == 1 ) { $gt = "0/1"; }
            elsif ( $maxGeno == 2 ) { $gt = "1/1"; }
            $$outfile_hash{$group}{gt} = $gt;

            my @filters = ();
            if ( $$outfile_hash{$group}{gq} < $opts{filter_quality} ) {
                push @filters, "q" . $opts{filter_quality};
            }
            if ( $numAlt < $opts{filter_evidence} ) {
                push @filters, "e" . $opts{filter_evidence};
            }
            if ( $numAlt + $numRef < $opts{filter_depth} ) {
                push @filters, "d" . $opts{filter_depth};
            }
            $$outfile_hash{$group}{ft} = ( defined( $filters[0] ) ) ? join( ";", @filters ) : "PASS";
        }
    }
    print "exiting scoreData()\n" if $opts{verbose};
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
    my ($outfile_hash) = @_;
    print "entering report()\n" if $opts{verbose};

    #open output file
    if ( defined( $opts{output_filename} ) ) {
        open( foutname1, ">$opts{output_filename}" ) or die("error opening file $opts{output_filename}\n");
    }
    else {
        open( foutname1, ">&", \*STDOUT ) or die;
    }

    if ( $opts{output_support} ) {

        #open support file
        open( support1, ">$opts{support_filename}" ) or die("could not open $opts{support_filename} for output, $!\n");
    }

    my $filedate = getDate();
    print foutname1 <<HEADER;
##fileformat=VCFv4.1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS:MT,Description="Nuclear Mitochondrial Insertion">
##FILTER=<ID=q$opts{filter_quality},Description="Phred-scaled quality filter">
##FILTER=<ID=e$opts{filter_evidence},Description="Support reads filter">
##FILTER=<ID=d$opts{filter_depth},Description="Sequence depth filter">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number likelihoods">
##FORMAT=<ID=CNL0,Number=.,Type=Float,Description="Copy number likelihoods with no frequency prior">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MSTART,Number=1,Type=Flag,Description="Mitochondrial start coordinate of inserted sequence">
##INFO=<ID=MEND,Number=1,Type=Flag,Description="Mitochondrial end coordinate of inserted sequence">
##INFO=<ID=MLEN,Number=1,Type=Flag,Description="Estimated length of mitochondrial insert">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##fileDate=$filedate
##reference=$opts{reference}
##source=dinumt-$version
HEADER

    my @vars = sort { $$outfile_hash{$a}{chr} cmp $$outfile_hash{$b}{chr} || $$outfile_hash{$a}{leftBkpt} <=> $$outfile_hash{$b}{leftBkpt} } keys %$outfile_hash;
    if ( $opts{output_gl} ) {
        print foutname1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$opts{prefix}\n";
    }
    else {
        print foutname1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    }

    my $index = 1;
    foreach my $group (@vars) {
        print "$group in report()\n" if $opts{verbose};
        if ( $$outfile_hash{$group}{gt} eq "0/0" || $$outfile_hash{$group}{gt} eq "./." ) {
            print "\thomref or nc, skipping\n" if $opts{verbose};
            next;
        }
        my $chrom = $$outfile_hash{$group}{chr};
        $chrom =~ s/chr//g;
        my $id      = $opts{prefix} . "_$index";
        my $alt     = "<INS:MT>";
        my $qual    = $$outfile_hash{$group}{gq};
        my $filter  = $$outfile_hash{$group}{ft};
        my %info    = ();
        my $pos     = $$outfile_hash{$group}{leftBkpt} - 1;
        my $end     = $$outfile_hash{$group}{rightBkpt};
        my $ciDelta = $end - $pos + 1;

        $info{IMPRECISE} = undef;
        $info{CIPOS}     = "0,$ciDelta";
        $info{CIEND}     = "-$ciDelta,0";
        $info{END}       = $end;
        $info{SVTYPE}    = "INS";

        if ( $$outfile_hash{$group}{m_len} ne "NA" ) {
            $info{MSTART} = $$outfile_hash{$group}{l_m_pos};
            $info{MEND}   = $$outfile_hash{$group}{r_m_pos};
            $info{MLEN}   = $$outfile_hash{$group}{m_len};
        }

        my $refline = `$opts{samtools} faidx $opts{reference} $chrom:$pos-$pos`;
        my $ref = ( split( /\n/, $refline ) )[1];
        if ( !defined($ref) ) { $ref = "N"; }

        my $format = "GT:FT:GL:GQ:PL";
        my $info   = "";
        my @sKeys  = sort { $a cmp $b } keys %info;
        for ( my $i = 0 ; $i <= $#sKeys ; $i++ ) {
            if ( $i > 0 ) { $info .= ";"; }
            if ( defined( $info{ $sKeys[$i] } ) ) {
                $info .= "$sKeys[$i]=$info{$sKeys[$i]}";
            }
            else {
                $info .= "$sKeys[$i]";
            }
        }
        if ( $opts{output_gl} ) {
            print foutname1 "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";
            my @gls = ();
            my @pls = ();
            foreach my $geno ( 0 .. $opts{ploidy} ) {
                push @gls, sprintf( "%.2f", $$outfile_hash{$group}{gl}{$geno} );
                push @pls, $$outfile_hash{$group}{pl}{$geno};
            }
            my $gl = join( ",", @gls );
            my $pl = join( ",", @pls );
            print foutname1 "\t$$outfile_hash{$group}{gt}:$$outfile_hash{$group}{ft}:$gl:$$outfile_hash{$group}{gq}:$pl\n";
        }
        else {
            print foutname1 "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\n";
        }
        if ( $opts{output_support} ) {
            print support1 "$$outfile_hash{$group}{support}";
        }
        $index++;
    }
    close(foutname1);
    if ( $opts{output_support} ) { close(support1); }

    print "exiting report()\n" if $opts{verbose};
}

sub calcGl {
    my ( $m, $g, $k, $l, $e ) = @_;
    print "in calcGl():\n"         if $opts{verbose};
    print "\t$m\t$g\t$k\t$l\t$e\n" if $opts{verbose};
    if ( 1 / $m**$k <= 0 ) { die "problem in calcGL 1, \t$m\t$g\t$k\t$l\t$e\n"; }
    my $gl = log10( 1 / $m**$k );
    if ( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) <= 0 ) { die "problem in calcGL 2, \t$m\t$g\t$k\t$l\t$e\n"; }
    $gl += log10( ( ( $m - $g ) * $e ) + ( ( 1 - $e ) * $g ) ) for 1 .. $l;
    if ( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) <= 0 ) { die "problem in calcGL 3, \t$m\t$g\t$k\t$l\t$e\n"; }
    $gl += log10( ( $m - $g ) * ( 1 - $e ) + ( $g * $e ) ) for ( $l + 1 ) .. $k;
    return $gl;
}

sub log10 {
    my $n = shift;
    return log($n) / log(10);
}

sub getInput {
    my ( $infile_hash, $readgroup_hash, $mask_hash, $mt_hash ) = @_;
    my @input_lines = ();
    print "Reading input files...\n" if $opts{verbose};

    #open input file
    if ( $opts{by_chr_dir} ) {
        if (defined( $opts{mt_names} ) ) { 
            foreach my $mt_name (keys %{$mt_hash}) {
                push @input_lines, "samtools view $opts{by_chr_dir}/$mt_name.*bam |";
            }
        }
        elsif ( $opts{ucsc} ) {
            push @input_lines, "samtools view $opts{by_chr_dir}/chrM.*bam |";
        }
        elsif ( $opts{ensembl} ) {
            push @input_lines, "samtools view $opts{by_chr_dir}/chrMT.*bam |";
        }
        else {
            push @input_lines, "samtools view $opts{by_chr_dir}/MT*bam |";
        }
    }
    else {
        if (defined( $opts{mt_names} ) ) {
            foreach my $mt_name (keys %{$mt_hash}) {
                push @input_lines, "samtools view $opts{input_filename} $mt_name |";
            }
        }
        elsif ( $opts{ucsc} ) {
            push @input_lines, "samtools view $opts{input_filename} chrM |";
        }
        elsif ( $opts{ensembl} ) {
            push @input_lines, "samtools view $opts{input_filename} chrMT |";
        }
        else {
            push @input_lines, "samtools view $opts{input_filename} MT |";
        }
    }

    #input mask coordinates
    if ( $opts{include_mask} ) {
        open( mask1, $opts{mask_filename} ) || die "Could not open $opts{mask_filename} for input, $!\n";
        while (<mask1>) {
            chomp;
            my ( $chr, $start, $end, $id ) = split(/\t/);
            $chr =~ s/chr//g;
            $$mask_hash{$chr}{$start} = $end;
            if ( $opts{include_mask} ) {
                if ( $opts{by_chr_dir} ) {
                    if ( $opts{ucsc} || $opts{ensembl} ) {
                        $chr = "chr" . $chr;        
                    }
                    push @input_lines, "samtools view $opts{by_chr_dir}/$chr*bam  $chr:$start-$end |";
                }
                else {
                    if ( $opts{ucsc} || $opts{ensembl} ) {  
                        $chr = "chr" . $chr;
                    }
                    push @input_lines, "samtools view $opts{input_filename} $chr:$start-$end |";
                }
            }
            $$mask_hash{$chr}{$start} = $end;
        }
        close mask1;
    }
    foreach my $input_line (@input_lines) {
        print "command: $input_line\n" if $opts{verbose};
        open( fname1, $input_line ) || die "error in opening file, $!\n";
        while ( my $line1 = <fname1> ) {
            $seq_num = $i;
            chomp($line1);
            my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split( /\t/, $line1 );

            my ($read_group) = $line1 =~ /RG:Z:(\S+)/;

            if ( $opts{read_groups} && !defined($read_group) ) { next; }
            elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

            if ( $rnext eq '=' || $rnext eq '*' ) { next; }
            if (defined( $opts{mt_names} ) ) {
                if (!defined($$mt_hash{$rname}) && defined($$mt_hash{$rnext}) ) { next; }
            }
            else {
                if ( $rname !~ /M/ && $rnext =~ /M/ ) { next; }
            }

            my $dnext = 0;
            if ( $flag & 32 ) {
                $dnext = 1;
            }

            my $dir = 0;
            if ( $flag & 16 ) {
                $dir = 1;
            }

            #Compare to masked regions
            my $isMaskOverlap = 0;
            if ( $opts{include_mask} ) {
                foreach my $maskStart ( keys %{ $$mask_hash{$rnext} } ) {
                    my $maskEnd = $$mask_hash{$rnext}{$maskStart};
                    if ( $pnext > $maskStart && $pnext <= $maskEnd ) {
                        $isMaskOverlap = 1;
                        last;
                    }
                }
            }
            if ($isMaskOverlap) { next; }

            #get mate information
            if ( $opts{by_chr_dir} ) {
                open( SAM, "samtools view $opts{by_chr_dir}/$rnext.*bam $rnext:$pnext-$pnext | " ) || die "Could not find MT bam file in $opts{by_chr_dir}, $!\n";
            }
            else {
                open( SAM, "samtools view $opts{input_filename} $rnext:$pnext-$pnext | " ) || die "Could not open $opts{input_filename}, $!\n";
            }

            my $c_mapq = 0;
            my $cnext  = 0;
            while (<SAM>) {
                chomp;
                my ( $m_qname, $m_flag, $m_rname, $m_pos, $m_mapq, $m_cigar, $m_rnext, $m_pnext, $m_tlen, $m_seq, $m_qual, $opt ) = split(/\t/);
                if ( $m_qname ne $qname ) { next; }
                $c_mapq = $m_mapq;
                $cnext  = $m_cigar;
                #if ( $c_mapq < $opts{min_map_qual} ) {
                #    print "MAPQ Filtering:\t$_\n" if $opts{verbose};
                #}
            }
            close SAM;

            if ( $c_mapq < $opts{min_map_qual} ) {
                next;
            }

            $$infile_hash{$seq_num}{group} = 0;
            $$infile_hash{$seq_num}{seq}   = $seq;

            $$infile_hash{$seq_num}{dir}   = $dir;
            $$infile_hash{$seq_num}{qname} = $qname;
            $$infile_hash{$seq_num}{rname} = $rname;
            $$infile_hash{$seq_num}{pos}   = $pos;
            $$infile_hash{$seq_num}{cigar} = $cigar;
            $$infile_hash{$seq_num}{cnext} = $cnext;
            $$infile_hash{$seq_num}{rnext} = $rnext;
            $$infile_hash{$seq_num}{pnext} = $pnext;
            $$infile_hash{$seq_num}{dnext} = $dnext;
            $$infile_hash{$seq_num}{tlen}  = $tlen;
            $$infile_hash{$seq_num}{qual}  = $qual;
            $$infile_hash{$seq_num}{line}  = $line1;

            $i++;
        }
        close(fname1);
    }
}

sub findBreakpoint {
    my ( $outfile_hash, $readgroup_hash, $mask_hash ) = @_;
    print "entering findBreakpoint()\n" if $opts{verbose};
    foreach my $group ( keys %$outfile_hash ) {
        print "\tAssessing group $group for breakpoints\n" if $opts{verbose};
        my $l_pos = $$outfile_hash{$group}{l_pos};
        my $r_pos = $$outfile_hash{$group}{r_pos};
        my $chro  = $$outfile_hash{$group}{chr};
        $$outfile_hash{$group}{numAltSR}  = 0;
        $$outfile_hash{$group}{numRefSR}  = 0;
        $$outfile_hash{$group}{leftBkpt}  = $l_pos;
        $$outfile_hash{$group}{rightBkpt} = $r_pos;

        #Compare to masked regions
        my $isMaskOverlap = 0;
        if ( $opts{include_mask} ) {
            foreach my $maskStart ( keys %{ $$mask_hash{$chro} } ) {
                my $maskEnd = $$mask_hash{$chro}{$maskStart};
                if ( ( $l_pos >= $maskStart && $l_pos <= $maskEnd ) || ( $r_pos >= $maskStart && $r_pos <= $maskEnd ) || ( $l_pos <= $maskStart && $r_pos >= $maskEnd ) ) {
                    $isMaskOverlap = 1;
                    last;
                }
            }
        }
        if ($isMaskOverlap) { next; }
        my %clippedPos = ();

        #open input file
        if ( $opts{by_chr_dir} ) {
            open( SAM, "samtools view $opts{by_chr_dir}/$chro.*bam $chro:$l_pos-$r_pos |" ) || die "Could not find MT bam file in $opts{by_chr_dir}, $!\n";
        }
        else {
            open( SAM, "samtools view $opts{input_filename} $chro:$l_pos-$r_pos |" ) || die "Could not open $opts{input_filename}, $!\n";
        }

        my %cnt = ();
        while (<SAM>) {
            chomp;
            my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);
            my ($read_group) = $_ =~ /RG:Z:(\S+)/;

            if ( $opts{read_groups} && !defined($read_group) ) { next; }
            elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

            #Check for read positions outside max_read_cov
            my $break = 0;
            for ( my $p = 0 ; $p <= length($seq) ; $p++ ) {
                $cnt{ $pos + $p }++;
                if ( $cnt{ $pos + $p } > $opts{max_read_cov} ) { $break = 1; last; }
            }

            if ($break) {
                print "Read count has reached limit of $opts{max_read_cov}, removing group $group\n" if $opts{verbose};
                delete( $$outfile_hash{$group} );
                last;
            }

            #Mark Clipped Positions
            my ( $cPos, $clipside, $clipsize ) = getSoftClipInfo( $pos, $cigar, $qual );
            if ( $cPos > -1 ) {
                $clippedPos{$cPos}++;
            }
        }
        close SAM;

        if ( !defined( $$outfile_hash{$group} ) ) { next; }
        my %bkpts    = ();
        my $numBkpts = 0;
        foreach my $cPos ( sort keys %clippedPos ) {
            if ( $clippedPos{$cPos} > 1 ) {
                $bkpts{$cPos} = $clippedPos{$cPos};
                $numBkpts++;
            }
        }

        my $num_bkpt_support = 0;
        my $num_ref_support  = 0;
        my $leftBkpt         = $l_pos;
        my $rightBkpt        = $r_pos;

        if ( $numBkpts > 0 && scalar keys %clippedPos <= $opts{max_num_clipped} ) {

            #take two most prevelant breaks for now
            my @sorted = sort { $bkpts{$b} <=> $bkpts{$a} } keys %bkpts;
            if ( $numBkpts == 1 ) {
                $leftBkpt         = $sorted[0];
                $rightBkpt        = $leftBkpt + 1;
                $num_bkpt_support = $bkpts{ $sorted[0] };
                $num_ref_support  = $cnt{ $sorted[0] } - $num_bkpt_support;
            }
            else {
                if ( $sorted[0] < $sorted[1] ) {
                    $leftBkpt  = $sorted[0];
                    $rightBkpt = $sorted[1];
                }
                else {
                    $leftBkpt  = $sorted[1];
                    $rightBkpt = $sorted[0];
                }
                $num_bkpt_support = $bkpts{ $sorted[0] } + $bkpts{ $sorted[1] };
                $num_ref_support  = $cnt{ $sorted[0] } + $cnt{ $sorted[1] } - $num_bkpt_support;
            }
        }
        $$outfile_hash{$group}{leftBkpt}  = $leftBkpt;
        $$outfile_hash{$group}{rightBkpt} = $rightBkpt;
        $$outfile_hash{$group}{numAltSR}  = $num_bkpt_support;
        $$outfile_hash{$group}{numRefSR}  = $num_ref_support;
    }
    print "exiting findBreakpoints()\n" if $opts{verbose};
}

sub seqCluster {
    my ($infile_hash) = @_;
    my $k             = 0;
    my %d             = ();

    $d{0}{k}     = 0;
    $d{0}{pnext} = 0;
    $d{0}{last}  = ();
    $d{0}{rnext} = 0;

    $d{1}{k}     = 0;
    $d{1}{pnext} = 0;
    $d{1}{last}  = ();
    $d{1}{rnext} = 0;

    my @sorted = sort { $$infile_hash{$a}->{rnext} cmp $$infile_hash{$b}->{rnext} || $$infile_hash{$a}->{pnext} <=> $$infile_hash{$b}->{pnext} } keys %{$infile_hash};
    print scalar @sorted . " total reads to process for clustering\n" if $opts{verbose};
    foreach my $c_seq_num (@sorted) {
        my $c_pnext = $$infile_hash{$c_seq_num}{pnext};
        my $c_dnext = $$infile_hash{$c_seq_num}{dnext};
        my $c_rnext = $$infile_hash{$c_seq_num}{rnext};
        my $c_qname = $$infile_hash{$c_seq_num}{qname};

        print "$c_dnext\n" if $opts{verbose};
        if ( $c_pnext - $d{$c_dnext}{pnext} > $opts{len_cluster_include} || $d{$c_dnext}{k} == 0 || $c_rnext ne $d{$c_dnext}{rnext} ) {

            #print "c_pnext:$c_pnext \t d_pnext:$d{$c_dnext}{pnext} \t dir:$d{$c_dnext}{k}\n" if $opts{verbose};
            if ( $d{$c_dnext}{k} > 0 && scalar @{ $d{$c_dnext}{last} } < $opts{min_reads_cluster} ) {
                foreach my $seq_num ( @{ $d{$c_dnext}{last} } ) {
                    delete $$infile_hash{$seq_num};
                }
            }
            $k++;
            $$infile_hash{$c_seq_num}{'group'} = $k;
            $d{$c_dnext}{k}                    = $k;
            $d{$c_dnext}{last}                 = ();
        }
        else {
            $$infile_hash{$c_seq_num}{'group'} = $d{$c_dnext}{k};

        }
        print "$d{$c_dnext}{k}\t$c_qname\t$c_rnext\t$c_pnext" if $opts{verbose};

        push @{ $d{$c_dnext}{last} }, $c_seq_num;
        $d{$c_dnext}{pnext} = $c_pnext;
        $d{$c_dnext}{rnext} = $c_rnext;

        #print "k:$k \t grp:$$infile_hash{$c_seq_num}{'group'} \t $d{$c_dnext}{k} \t pre_pos:$d{$c_dnext}{pnext}\n" if $opts{verbose};
        #print "$$infile_hash{$c_seq_num}->{'group'}\n chr_num:$c_rnext\n"                                          if $opts{verbose};
    }

}

sub linkCluster {
    my ($infile_hash) = @_;

    #this can link multiple F's to a single leftmost R
    my @sorted = sort { $$infile_hash{$a}->{rnext} cmp $$infile_hash{$b}->{rnext} || $$infile_hash{$a}->{pnext} <=> $$infile_hash{$b}->{pnext} } keys %{$infile_hash};
    print scalar @sorted . " total reads to process for linking clusters\n" if $opts{verbose};

    for ( my $c = 0 ; $c <= $#sorted ; $c++ ) {
        my $c_seq_num = $sorted[$c];
        $$infile_hash{$c_seq_num}{link} = 0;
        my $c_pnext = $$infile_hash{$c_seq_num}{pnext};
        my $c_dnext = $$infile_hash{$c_seq_num}{dnext};
        my $c_rnext = $$infile_hash{$c_seq_num}{rnext};
        my $c_dir   = $$infile_hash{$c_seq_num}{dir};

        if ( $c_dnext == 1 ) { next; }

        for ( my $d = $c + 1 ; $d <= $#sorted ; $d++ ) {
            my $d_seq_num = $sorted[$d];
            my $d_pnext   = $$infile_hash{$d_seq_num}{pnext};
            my $d_dnext   = $$infile_hash{$d_seq_num}{dnext};
            my $d_rnext   = $$infile_hash{$d_seq_num}{rnext};
            my $d_dir     = $$infile_hash{$d_seq_num}{dir};

            if ( $d_dnext == 0 ) { next; }

            my $delta = $d_pnext - $c_pnext;

            if (
                   $delta < $opts{len_cluster_link}
                && $c_rnext eq $d_rnext
                && (   ( $c_dir == 0 && $c_dnext == 1 && $d_dnext == 0 && $d_dir == 1 )
                    || ( $c_dir == 0 && $c_dnext == 0 && $d_dnext == 1 && $d_dir == 1 )
                    || ( $c_dir == 1 && $c_dnext == 0 && $d_dnext == 1 && $d_dir == 0 ) )
              )
            {

                #0.0.11 - must have consistent orientation between left and right sides of insertions
                $$infile_hash{$c_seq_num}{link} = $$infile_hash{$d_seq_num}{group};
                $$infile_hash{$d_seq_num}{link} = $$infile_hash{$c_seq_num}{group};
            }
            elsif ( !defined( $$infile_hash{$c_seq_num}{link} ) ) {
                $$infile_hash{$c_seq_num}{link} = 0;
            }
        }
    }
}

sub mapCluster {
    my ( $infile_hash, $outfile_hash, $readgroup_hash ) = @_;

    my %l_linked_groups;
    my %linked_group_pnext;

    foreach my $key ( sort { $infile_hash{$a}->{group} <=> $infile_hash{$b}->{group} } keys %infile_hash ) {
        if ( ( $infile_hash{$key}{'group'} > 0 ) && ( $infile_hash{$key}{'dnext'} == 0 ) && ( $infile_hash{$key}{'link'} > 0 ) ) {
            $l_linked_groups{ $infile_hash{$key}{'group'} } = $infile_hash{$key}{'link'};
        }
    }

    my $i = 1;
    while ( my ( $group, $link ) = each(%l_linked_groups) ) {
        print "group = $group ;; link = $link \n" if $opts{verbose};
        my @l_rnext = map { $infile_hash{$_}{'rnext'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_rname = map { $infile_hash{$_}{'rname'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_pos   = map { $infile_hash{$_}{'pos'} } grep   { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_dir   = map { $infile_hash{$_}{'dir'} } grep   { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_qname = map { $infile_hash{$_}{'qname'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_pnext = map { $infile_hash{$_}{'pnext'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_dnext = map { $infile_hash{$_}{'dnext'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_cigar = map { $infile_hash{$_}{'cigar'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @l_cnext = map { $infile_hash{$_}{'cnext'} } grep { $infile_hash{$_}{group} == $group } keys %infile_hash;
        my @r_rnext = map { $infile_hash{$_}{'rnext'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_pos   = map { $infile_hash{$_}{'pos'} } grep   { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_dir   = map { $infile_hash{$_}{'dir'} } grep   { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_rname = map { $infile_hash{$_}{'rname'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_pnext = map { $infile_hash{$_}{'pnext'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_dnext = map { $infile_hash{$_}{'dnext'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_qname = map { $infile_hash{$_}{'qname'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_cigar = map { $infile_hash{$_}{'cigar'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;
        my @r_cnext = map { $infile_hash{$_}{'cnext'} } grep { $infile_hash{$_}{group} == $link } keys %infile_hash;

        my @rc_pnext = ();
        my @lc_pnext = ();
        my $chr      = $l_rnext[0];

        #update right coordinates based on cigar length
        for ( my $c = 0 ; $c <= $#r_cnext ; $c++ ) {
            my $cigar = $r_cnext[$c];
            $rc_pnext[$c] = $r_pnext[$c];

            while ( $cigar =~ /(\d+)M/g ) {
                $rc_pnext[$c] += $1;
            }
            while ( $cigar =~ /(\d+)N/g ) {
                $rc_pnext[$c] += $1;
            }
            while ( $cigar =~ /(\d+)D/g ) {
                $rc_pnext[$c] += $1;
            }
        }

        #update left coordinates based on cigar length
        for ( my $c = 0 ; $c <= $#l_cnext ; $c++ ) {
            my $cigar = $l_cnext[$c];
            $lc_pnext[$c] = $l_pnext[$c];

            while ( $cigar =~ /(\d+)M/g ) {
                $lc_pnext[$c] += $1;
            }
            while ( $cigar =~ /(\d+)N/g ) {
                $lc_pnext[$c] += $1;
            }
            while ( $cigar =~ /(\d+)D/g ) {
                $lc_pnext[$c] += $1;
            }
        }

        my @s_l_pnext = sort { $a <=> $b } @l_pnext;
        my @s_r_pnext = sort { $a <=> $b } @r_pnext;

        my @s_lc_pnext = sort { $a <=> $b } @lc_pnext;
        my @s_rc_pnext = sort { $a <=> $b } @rc_pnext;

        my $l_brk_point = $s_l_pnext[$#s_l_pnext];
        my $r_brk_point = $s_rc_pnext[0];

        my $win_l_s = $s_l_pnext[0];
        my $win_l_e = $s_lc_pnext[$#s_lc_pnext];
        my $win_r_s = $s_r_pnext[0];
        my $win_r_e = $s_rc_pnext[$#s_rc_pnext];

        $linked_group_pnext{$i}{l_group} = $group;
        $linked_group_pnext{$i}{r_group} = $link;

        #Check for crossing clusters; occurs in regions of high coverage, but not ascertained
        #until next step so stopgap here
        if ( $l_brk_point > $r_brk_point ) {
            my $temp = $l_brk_point;
            $l_brk_point = $r_brk_point;
            $r_brk_point = $temp;
        }

        my $commandLeft  = "";
        my $commandRight = "";
        if ( $opts{by_chr_dir} ) {
            $commandLeft  = "samtools view $opts{by_chr_dir}/$chr.*bam $chr:$win_l_s-$win_l_e |";
            $commandRight = "samtools view $opts{by_chr_dir}/$chr.*bam $chr:$win_r_s-$win_r_e |";
        }
        else {
            $commandLeft  = "samtools view $opts{input_filename} $chr:$win_l_s-$win_l_e |";
            $commandRight = "samtools view $opts{input_filename} $chr:$win_r_s-$win_r_e |";
        }

        my $numRefRP = 0;
        my $numAltRP = scalar @l_pnext + scalar @r_pnext;
        my $sumE     = 0;

        #left
        open( SAM, $commandLeft ) || die "Could not open sam file for input, $!\n";
        while (<SAM>) {
            chomp;
            my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);

            my $mapE = 10**( -1 * $mapq / 10 );
            my ($read_group) = $_ =~ /RG:Z:(\S+)/;

            if ( $opts{read_groups} && !defined($read_group) ) { next; }
            elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

            if ( $mapq < $opts{min_map_qual} ) { next; }
            my $dir = 0;    #F
            if ( $flag & 16 ) { $dir = 1; }    #R

            if ( $pos >= $win_l_s && $pos <= $win_l_e && $dir == 0 ) {
                $numRefRP++;
                $sumE += $mapE;
            }
        }
        close SAM;

        #right
        open( SAM, $commandRight ) || die "Could not open sam file for input, $!\n";
        while (<SAM>) {
            chomp;
            my ( $qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split(/\t/);

            my $mapE = 10**( -1 * $mapq / 10 );
            my ($read_group) = $_ =~ /RG:Z:(\S+)/;

            if ( $opts{read_groups} && !defined($read_group) ) { next; }
            elsif ( $opts{read_groups} && !defined( $$readgroup_hash{$read_group} ) ) { next; }

            if ( $mapq < $opts{min_map_qual} ) { next; }
            my $dir = 0;    #F
            if ( $flag & 16 ) { $dir = 1; }    #R

            if ( $pos >= $win_r_s && $pos <= $win_r_e && $dir == 1 ) {
                $numRefRP++;
                $sumE += $mapE;
            }
        }
        close SAM;

        my $avgQ = $sumE / $numRefRP;          #calculate average over all reads, ref and alt
        $numRefRP -= $numAltRP;                #correct for alt reads

        #estimate mitochondrial coordinates from mated sequence alignments
        my $l_m_min   = 1e10;
        my $l_m_max   = 0;
        my $l_m_min_i = -1;
        my $l_m_max_i = -1;
        my $l_n_dir   = -1;
        my $l_m_dir   = -1;

        for ( my $i = 0 ; $i <= $#l_qname ; $i++ ) {
            if ( $l_rname[$i] !~ /M/ ) { next; }    #don't include nuclear homologous regions
            if ( $l_pos[$i] < $l_m_min ) {
                $l_m_min   = $l_pos[$i];
                $l_m_min_i = $i;
            }
            if ( $l_pos[$i] > $l_m_max ) {
                $l_m_max   = $l_pos[$i];
                $l_m_max_i = $i;
            }
            $l_n_dir = $l_dnext[$i];
            $l_m_dir = $l_dir[$i];
        }
        my $r_m_min   = 1e10;
        my $r_m_max   = 0;
        my $r_m_min_i = -1;
        my $r_m_max_i = -1;
        my $r_n_dir   = -1;
        my $r_m_dir   = -1;

        for ( my $i = 0 ; $i <= $#r_qname ; $i++ ) {
            if ( $r_rname[$i] !~ /M/ ) { next; }    #don't include nuclear homologous regions
            if ( $r_pos[$i] < $r_m_min ) {
                $r_m_min   = $r_pos[$i];
                $r_m_min_i = $i;
            }
            if ( $r_pos[$i] > $r_m_max ) {
                $r_m_max   = $r_pos[$i];
                $r_m_max_i = $i;
            }
            $r_n_dir = $r_dnext[$i];
            $r_m_dir = $r_dir[$i];
        }
        $$outfile_hash{$group}{l_m_pos} = "NA";
        $$outfile_hash{$group}{r_m_pos} = "NA";
        $$outfile_hash{$group}{m_len}   = "NA";

        if ( $l_m_dir > -1 && $r_m_dir > -1 ) {    #have mitochondrial mappings
            if ( $l_n_dir == 0 && $l_m_dir == 1 && $r_m_dir == 0 && $r_n_dir == 1 ) {
                my $cigar = $r_cigar[$r_m_max_i];
                while ( $cigar =~ /(\d+)M/g ) {
                    $r_m_max += $1;
                }
                while ( $cigar =~ /(\d+)N/g ) {
                    $r_m_max += $1;
                }
                while ( $cigar =~ /(\d+)D/g ) {
                    $r_m_max += $1;
                }
            }
            elsif ( $l_n_dir == 0 && $l_m_dir == 0 && $r_m_dir == 1 && $r_n_dir == 1 ) {
                my $cigar = $l_cigar[$l_m_max_i];
                while ( $cigar =~ /(\d+)M/g ) {
                    $l_m_max += $1;
                }
                while ( $cigar =~ /(\d+)N/g ) {
                    $l_m_max += $1;
                }
                while ( $cigar =~ /(\d+)D/g ) {
                    $l_m_max += $1;
                }
            }

            $$outfile_hash{$group}{l_m_pos} = $l_m_min;
            $$outfile_hash{$group}{r_m_pos} = $r_m_max;
            if ( $r_m_max > $l_m_min ) {
                $$outfile_hash{$group}{l_m_pos} = $l_m_min;
                $$outfile_hash{$group}{r_m_pos} = $r_m_max;
            }
            else {
                $$outfile_hash{$group}{l_m_pos} = $r_m_max;
                $$outfile_hash{$group}{r_m_pos} = $l_m_max;
            }

            #currently assumes smallest sequence possible due to circular nature of mt dna
            $$outfile_hash{$group}{m_len} = $$outfile_hash{$group}{r_m_pos} - $$outfile_hash{$group}{l_m_pos} + 1;
            my $lenAlt = $opts{len_mt} - $$outfile_hash{$group}{r_m_pos} + $$outfile_hash{$group}{l_m_pos} + 1;
            if ( $lenAlt < $$outfile_hash{$group}{m_len} ) { $$outfile_hash{$group}{m_len} = $lenAlt; }
        }

        $$outfile_hash{$group}{avgQ}     = $avgQ;
        $$outfile_hash{$group}{l_pos}    = $l_brk_point;
        $$outfile_hash{$group}{r_pos}    = $r_brk_point;
        $$outfile_hash{$group}{chr}      = $chr;
        $$outfile_hash{$group}{numRefRP} = $numRefRP;
        $$outfile_hash{$group}{numAltRP} = $numAltRP;

        if ( $opts{output_support} ) {
            foreach my $key ( grep { $infile_hash{$_}{group} == $group } keys %infile_hash ) {
                $$outfile_hash{$group}{support} .= "$infile_hash{$key}{line}\n";
            }
            foreach my $key ( grep { $infile_hash{$_}{group} == $link } keys %infile_hash ) {
                $$outfile_hash{$group}{support} .= "$infile_hash{$key}{line}\n";
            }
        }
        $i++;

        if ( $opts{verbose} ) {
            print "\t$chr\t$$outfile_hash{$group}{l_pos}\t$$outfile_hash{$group}{r_pos}\t$$outfile_hash{$group}{numRefRP}\t$$outfile_hash{$group}{numAltRP}\t$$outfile_hash{$group}{avgQ}\t$$outfile_hash{$group}{l_m_pos}\t$$outfile_hash{$group}{r_m_pos}\t$$outfile_hash{$group}{m_len}\n";
        }
    }
}

sub getSoftClipInfo {
    my ( $pos, $cigar, $qual ) = @_;
    my $clipside = "";
    my $clipsize = 0;
    my $cPos     = -1;
    my $avgQual  = -1;

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
            $clippedQuals = substr( $qual, length($qual) - $clipsize - 1, $clipsize );
        }
        else {
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
    if ( $avgQual < 10 ) { $cPos = -1; }
    if ( $clipsize < $opts{min_clipped_seq} ) { $cPos = -1; }
    return ( $cPos, $clipside, $clipsize );
}

sub usage {
    my $version = shift;
    printf("\n");
    printf( "%-9s %s\n", "Program:", "dinumt.pl" );
    printf( "%-9s %s\n", "Version:", "$version" );
    printf("\n");
    printf( "%-9s %s\n", "Usage:", "dinumt.pl [options]" );
    printf("\n");
    printf( "%-9s %-35s %s\n", "Options:", "--input_filename=[filename]",     "Input alignment file in BAM format" );
    printf( "%-9s %-35s %s\n", "",         "--output_filename=[filename]",    "Output file (default stdout)" );
    printf( "%-9s %-35s %s\n", "",         "--mask_filename=[filename]",      "Mask file for reference numts in BED format (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--include_mask",                  "Include aberrant reads mapped to mask regions in clustering" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_include=[integer]", "Maximum distance to be included in cluster (default 600)" );
    printf( "%-9s %-35s %s\n", "",         "--len_cluster_link=[integer]",    "Maximum distance to link clusters (default 800)" );
    printf( "%-9s %-35s %s\n", "",         "--min_reads_cluster=[integer]",   "Minimum number of reads to link a cluster (default 1)" );
    printf( "%-9s %-35s %s\n", "",         "--min_evidence=[integer]",        "Minimum evidence to consider an insertion event (default 4)" );
    printf( "%-9s %-35s %s\n", "",         "--min_map_qual=[integer]",        "Minimum mapping quality for read consideration (default 10)" );
    printf( "%-9s %-35s %s\n", "",         "--max_read_cov=[integer]",        "Maximum read coverage allowed for breakpoint searching (default 200)" );
    printf( "%-9s %-35s %s\n", "",         "--min_clipped_seq=[integer]",     "Minimum clipped sequence required to consider as putative breakpoint (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--max_num_clipped=[integer]",     "Maximum number of clipped sequences observed before removing from evidence consideration (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--read_groups=[read_group1],...", "Limit analysis to specified read group(s)" );
    printf( "%-9s %-35s %s\n", "",         "--mt_names=[mt_name1],...",       "Limit analysis to specified mitochondrial sequence names" );
    printf( "%-9s %-35s %s\n", "",         "--by_chr_dir=[directory]",        "If set, expects to find chr specific BAM files in indicated directory" );
    printf( "%-9s %-35s %s\n", "",         "--prefix=[string]",               "Prepend label in report output" );
    printf( "%-9s %-35s %s\n", "",         "--ucsc",                          "Use UCSC genome formatting (e.g. chrM)" );
    printf( "%-9s %-35s %s\n", "",         "--ensembl",                       "Use Ensembl genome formatting (e.g. chrMT)" );
    printf( "%-9s %-35s %s\n", "",         "--output_gl",                     "Output genotype likelihood information" );
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
    if ( $$opts{include_mask} && !defined( $$opts{mask_filename} ) ) {
        print "\n***ERROR***\t--mask_filename is neccessary with --include_mask option\n";
        usage($version);
        exit;
    }
    if ( $$opts{output_support} && !defined( $$opts{output_filename} ) ) {
        print "\n***ERROR***\t--output_filename is neccessary with --output_support option\n";
        usage($version);
        exit;
    }
    if ( $$opts{ucsc} && $$opts{ensembl} ) {
        print "\n***ERROR***\t--ucsc and --ensembl are mutually exclusive options\n";
        usage($version);
        exit;
    }
}
