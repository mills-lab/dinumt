#!/usr/bin/perl

use strict;
use warnings;

my %data = ();

my $lastChr = "";
my $lastPos = 0;
my $index = 0;
my $cnt = 1;
my $last = "";
my $sumLen = 0;
my $sumSupport =0;
my $c = 0;
my $len_mt = 16596;

while (<>) {
    chomp;
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split(/\t/);
    if (($chr eq $lastChr && $pos - $lastPos > 1000) || ($chr ne $lastChr && $lastChr ne "")) {
        $c++;
        $index = 0;
    }
    my ($end) = $info =~ /END=(\d+)/;
    my ($mlen) = $info =~ /MLEN=(\d+)/;
    my ($mstart) = $info =~ /MSTART=(\d+)/;
    my ($mend) = $info =~ /MEND=(\d+)/;
    $data{$c}[$index]{chr} = $chr;
    $data{$c}[$index]{pos} = $pos;
    $data{$c}[$index]{end} = $end;
    $data{$c}[$index]{id} = $id;
    $data{$c}[$index]{qual} = $qual;
    $data{$c}[$index]{filter} = $filter;
    $data{$c}[$index]{mlen} = (defined($mlen)) ? $mlen : "NA";
    $data{$c}[$index]{mstart} = (defined($mstart)) ? $mstart : "NA";
    $data{$c}[$index]{mend} = (defined($mend)) ? $mend : "NA";

    $lastChr = $chr;

    $lastPos = $pos;
    $index++;
}

my $mergedId = 0;
for (my $i=0; $i<=$c; $i++) {
    my $n = 0;
    my %seg = ();
    
    my @sortedStarts = sort { $a->{pos} <=> $b->{pos} } @{$data{$i}};
    my @sortedEnds = sort { $a->{end} <=> $b->{end} } @{$data{$i}};

    if ($#sortedStarts == 0 && $data{$i}[0]{id} =~ /DG196/) { next; } 
    
    my $qual_sum = 0;
    my $mstart_min = 1e100;
    my $mend_max = 0;
    my %samples = ();
    my $filter = "LowQual";
    my $chr = $sortedStarts[0]->{chr};

    for (my $s=0; $s<=$#sortedStarts; $s++) {
        my $sample = $sortedStarts[$s]->{id};
        if ($sortedStarts[$s]->{mstart} ne "NA") {
            if ($sortedStarts[$s]->{mstart} < $mstart_min) { $mstart_min = $sortedStarts[$s]->{mstart}; }
            if ($sortedStarts[$s]->{mend} > $mend_max) { $mend_max = $sortedStarts[$s]->{mend}; }
        }

        $sample =~ s/_.*//;
        $samples{$sample} = 1;
        $qual_sum += $sortedStarts[$s]->{qual};
        if ($sortedStarts[$s]->{filter} eq "PASS") { $filter = "PASS"; }  
        my $lastD = -1;
        if ($s > 0) { $lastD = $sortedStarts[$s - 1]->{pos}; } 
        for (my $e=0; $e < $s; $e++) {
            if ($sortedEnds[$e]->{end} >= $sortedStarts[$s - 1]->{pos} && 
                $sortedEnds[$e]->{end} <= $sortedStarts[$s]->{pos} ) {
                    if (!defined($seg{$lastD})) {
                        $seg{$lastD}{end} = $sortedEnds[$e]->{end};
                        $seg{$lastD}{n} = $n;
                    }
                    $n--;
                    $lastD = $sortedEnds[$e]->{end};
            }
        }
        if ($lastD > -1) {
            $seg{$lastD}{end} = $sortedStarts[$s]->{pos};
            $seg{$lastD}{n} = $n;
        }
        $n++; 
    }

    if (scalar @sortedStarts == 1) { $seg{$sortedStarts[0]->{pos}}{end} = $sortedStarts[0]->{end};  $seg{$sortedStarts[0]->{pos}}{n} = 1; } 
    $mergedId++;
    my %info = ();

    my @sortedSeg = sort { $seg{$b}{n} <=> $seg{$a}{n} } keys %seg; 

    my $id = "MERGED_NUMT_$mergedId";
    my $total = scalar @sortedStarts;
    my $qual = int($qual_sum / $total);
    my $pos = $sortedSeg[0];
    my $end = $seg{$sortedSeg[0]}{end};
    my $ciDelta = $end - $pos + 1;
    my $alt     = "<INS:MT>";

    #my $refline = `samtools faidx /mnt/Mills/reference/GRCh37/Sequence/WholeGenomeFasta/genome.fa $chr:$pos-$pos`;
    my $refline = `samtools faidx /mnt/Mills/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa $chr:$pos-$pos`;
    my $ref = ( split( /\n/, $refline ) )[1];
    if ( !defined($ref) ) { $ref = "N"; }

    $info{END} = $end; 
    $info{SAMPLES} = join(",", sort keys %samples); 
    $info{IMPRECISE} = undef;
    $info{CIPOS}     = "0,$ciDelta";
    $info{CIEND}     = "-$ciDelta,0";
    $info{SVTYPE}    = "INS";
    if ($mend_max > 0) {
        $info{MSTART} = $mstart_min;
        $info{MEND} = $mend_max;
        $info{MLEN} = ($mend_max - $mstart_min + 1);
        if (abs($mend_max - $len_mt + $mstart_min  + 1) < $info{MLEN}) { $info{MLEN} = abs($mend_max - $len_mt + $mstart_min + 1); } 
    }
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
    
    print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\n";
}
