#!/usr/bin/perl -w
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#                  PRC, the profile comparer, version 1.5.6
#
#  merge_aligns.pl: Perl script for combining HMM-seq and HMM-HMM alignments
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#  Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
#  All Rights Reserved
#
#  This source code is distributed under the terms of the GNU General Public 
#  License. See the files COPYING and LICENSE for details.
#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
use strict;

my $about = "
Usage:

  merge_aligns.pl hmm1-hmm2.out hmm1-seq1.a2m hmm2-seq2.a2m > seq1-seq2.a2m

This script merges an HMM-HMM alignment created by PRC using 

  prc -hits 1 -align prc hmm1.mod hmm2.mod > hmm1-hmm2.out 

with SAM HMM-seq alignemnts created using e.g.

  hmmscore hmm1-seq1 -i hmm1.mod -db seq1.fa \\
    -sw 2 -dpstyle 0 -adpstyle 5 -select_align 8

  hmmscore hmm2-seq2 -i hmm2.mod -db seq2.fa \\
    -sw 2 -dpstyle 0 -adpstyle 5 -select_align 8

The resulting alignment is reported in the .a2m format (aligned
residues are uppercase, unaligned residues are lowercase).

This script uses the first alignment in the PRC output file, which is
also the one with the highest reverse score.

";

die $about unless scalar(@ARGV)==3;


# parse a FASTA file (also works for .a2m files)
# returns something like:
#
#    [ { seqid => 'seqid1',
#        seq   => 'ASFDSASDFASDFASDFASDF' }, 
#
#      { seqid => 'seqid2',
#        seq   => 'ASFASFASDFASDFASDF' } ]
#
sub read_fasta($)
{
    my ($file) = @_;
    my (@seqs, $r);

    open F, $file
        or die "Cannot open '$file': $!\n ... ";

  Main:
    while(<F>)
    {
        next unless /^>\s*(\S+)/;

	$r = { seqid => $1, seq => '' };
	push @seqs, $r;

        while(<F>)
        {
            redo Main if /^>/;
            next unless /^\s*(\S+)/;
            $r->{seq} .= $1;
        }
    }
    close F;

    return \@seqs;
}

# NB this routine also works for .a2m sequences
sub print_fasta($$)
{
    my ($seqid, $seq) = @_;
    
    print ">$seqid\n";

    while(length($seq)>50)
    {
	$seq =~ s/^(\S{50})//;
	print "$1\n";
    }

    print "$seq\n";
}

# read a PRC alignment generated using
#
#   $ prc -hits 1 -align prc hmm1 hmm2 > file.out
#
sub read_prc($)
{
    my ($prc_file) = @_;
    my ($prc_al1, $prc_al2);

    # get the two alignment strings
    open F, $prc_file
	or die "Cannot open '$prc_file' for reading: $!\n ... ";

    while(<F>){ last if /^\#\shmm1/ };
    my ($hmm1_id, $start1, $end1, $length1, $hit_no,
	$hmm2_id, $start2, $end2, $length2) = split /\t/, <F>;

    die "Error parsing the hit of $prc_file!\n ... "
	unless defined $length2;

    die "The first hit_no isn't 1!\n ... "
	unless $hit_no==1;

    <F>;
    <F> =~ /^([MDI~]+)/
	or die "Error parsing alignment line 1 of $prc_file!\n ... ";
    $prc_al1 = $1;

    <F> =~ /^([MDI~]+)/
	or die "Error parsing alignment line 2 of $prc_file!\n ... ";
    $prc_al2 = $1;

    close F;

    die "The two alignment lines are of different lengths!\n ... "
	unless length($prc_al1) == length($prc_al2);

    # prepend and append xxx for unmatched parts of the model
    $prc_al1 = 'x' x ($start1-1) . $prc_al1 . 'x' x ($length1-$end1);
    $prc_al2 = 'x' x ($start2-1) . $prc_al2 . 'x' x ($length2-$end2);

    return ( $prc_al1, $prc_al2 );
}

# merge HMM1-HMM2 alignments in the PRC two-line format with
# HMM1-seed1 and HMM2-seed2 alignments (in the usual SAM format), and
# print out the seed1-seed2 a2m alignment
#
sub merge_prc_seed($$$$$$)
{
    my ($chain1_id, $chain2_id, $prc_al1, $prc_al2, $seed_al1, $seed_al2) = @_;

    # first check that the lengths are sane
    {
	my ($p1, $p2, $s1, $s2) = ($prc_al1, $prc_al2, $seed_al1, $seed_al2);

	my $lp1 = ($p1 =~ tr/[xMD]//);
	my $lp2 = ($p2 =~ tr/[xMD]//);
	my $ls1 = ($s1 =~ tr/[A-Z\-]//);
	my $ls2 = ($s2 =~ tr/[A-Z\-]//);

	die "The PRC/SAM lengths don't match for chain1 ($lp1/$ls1)!\n ... "
	    unless $lp1 == $ls1;

	die "The PRC/SAM lengths don't match for chain2 ($lp2/$ls2)!\n ... "
	    unless $lp2 == $ls2;
    }


    my @seed1 = split //, $seed_al1;
    my @seed2 = split //, $seed_al2;
    my @prc1 = split //, $prc_al1;
    my @prc2 = split //, $prc_al2;

    my ($a2m1, $a2m2, $seed1, $seed2, $prc1, $prc2);
    $a2m1 = ''; $a2m2 = '';

  Main:
    while( scalar(@seed1)>0 or scalar(@seed2)>0 )
    {
	#
	# first make sure $seed1,2 and $prc1,2 all defined
	#

	while( scalar(@seed1)>0 and not defined $seed1 )
	{
	    $seed1 = shift @seed1;

	    if( $seed1 =~ /^[a-z]$/ )
	    {
		# unaligned seed residues will always remain unaligned
		# so skip right through them
		$a2m1 .= $seed1;
		undef $seed1;
		next Main;
	    }   
	    elsif( $seed1 !~ /^[A-Z\-]$/ )
	    {
		# skip all weird characters
		undef $seed1;
		next Main;
	    }
	}

	while( scalar(@seed2)>0 and not defined $seed2 )
	{
	    $seed2 = shift @seed2;

	    if( $seed2 =~ /^[a-z]$/ )
	    {
		$a2m2 .= $seed2;
		undef $seed2;
		next Main;
	    }   
	    elsif( $seed2 !~ /^[A-Z\-]$/ )
	    {
		undef $seed2;
		next Main;
	    }
	}
	
	$prc1 = shift @prc1 if scalar(@prc1)>0 and not defined $prc1;
	$prc2 = shift @prc2 if scalar(@prc2)>0 and not defined $prc2;

	# 'x' = region of the model outside of PRC match
	# skip right through

	if(defined $prc1 and $prc1 eq 'x')
	{
	    die "Uh-oh, chain1 alignment lengths don't match up!\n ... "
		unless defined $seed1;

	    $a2m1 .= lc $seed1 if $seed1 ne '-';
	    undef $seed1; undef $prc1;
	    next Main;
	}

	if(defined $prc2 and $prc2 eq 'x')
	{
	    die "Uh-oh, chain2 alignment lengths don't match up!\n ... "
		unless defined $seed2;

	    $a2m2 .= lc $seed2 if $seed2 ne '-';
	    undef $seed2; undef $prc2;
	    next Main;
	}

	#print "$prc1($seed1) -- $prc2($seed2)\n";
	
	die "Something's screwed up!\n ... "
	    if( not defined $prc1 or not defined $prc2 
		or ($prc1 ne '~' and $prc1 ne 'I' and not defined $seed1)
		or ($prc2 ne '~' and $prc2 ne 'I' and not defined $seed2) );

	#
	# now figure out what to do about the model-model alignment
	#

	if( $prc1 eq 'M' and $prc2 eq 'M' )
	{
	    if( $seed1 ne '-' and $seed2 ne '-' )
	    {
		# proper match!!!!
		$a2m1 .= uc $seed1;
		$a2m2 .= uc $seed2;
	    }
	    elsif( $seed1 ne '-' )
	    {
		# $seed2 is '-', so nothing to align $seed1 to!
		$a2m1 .= lc $seed1;
	    }
	    elsif( $seed2 ne '-' )
	    {
		# $seed1 is '-', so nothing to align $seed2 to!
		$a2m2 .= lc $seed2;
	    }
	    else
	    {
		# both are '-': do nothing!
		;
	    }

	    undef $seed1; undef $prc1;
	    undef $seed2; undef $prc2;
	}
	elsif( ( $prc1 eq 'M' and $prc2 eq 'I' ) or
	       ( $prc1 eq 'D' and $prc2 eq '~' )    )
	{
	    $a2m1 .= lc $seed1 if $seed1 ne '-';

	    undef $seed1; undef $prc1;
	    undef $prc2;
	}
	elsif( ( $prc1 eq 'I' and $prc2 eq 'M' ) or
	       ( $prc1 eq '~' and $prc2 eq 'D' )    )
	{
	    $a2m2 .= lc $seed2 if $seed2 ne '-';

	    undef $prc1;
	    undef $seed2; undef $prc2;
	}
	else
	{
	    die "Wrong logic somewhere!\n ... ";
	}	
    }

    print_fasta($chain1_id, $a2m1);
    print_fasta($chain2_id, $a2m2);
}

# main
{
    my ($prc_file, $seed1_file, $seed2_file) = @ARGV;

    my $file1 = read_fasta($seed1_file);
    my ($chain1_id, $seed1_al) = ( $file1->[0]{seqid}, $file1->[0]{seq} );
    
    my $file2 = read_fasta($seed2_file);
    my ($chain2_id, $seed2_al) = ( $file2->[0]{seqid}, $file2->[0]{seq} );
    
    my ($prc1_al, $prc2_al) 
	= read_prc($prc_file);
    
    merge_prc_seed( $chain1_id, $chain2_id, 
		    $prc1_al, $prc2_al, $seed1_al, $seed2_al );
}
