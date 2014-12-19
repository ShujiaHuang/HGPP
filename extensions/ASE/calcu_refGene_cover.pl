# Author : Shujia Huang
# Date   : 2012-07-09
#!/usr/bin/perl
use strict;
use warnings;

my ( $geneRegion, $scaffoldId, $bamfileT, $bamfileHomo, $bamfileAmb, $bamfileOther ) = @ARGV;

die "perl $0 [gene_region] [scaffoldId] [bamfileT] <[bamfileHomo] [bamfileAmb] [bamfileOther]>\n" if ( @ARGV < 3 );

my ( %cds, %gene, %cvgReg, %index );

my @tmp;
open I, $geneRegion or die "Cannot open file : $geneRegion\n";
while ( <I> ) {
#scaffold238     mRNA    156154  179047  -       NM_012346       NUP62   chr19,mRNA,50410083,50432773,-,NM_012346,NUP62
#scaffold238     CDS     176030  177599  -       NM_012346       NUP62   chr19,CDS,50411495,50413064,-,NM_012346,NUP62
	chomp;
	@tmp    = split;
	my $key = $tmp[0].":".$tmp[6];
	next if ( $tmp[1] ne "mRNA" ); # Just use mRNA
	# Caution: there's not all the genes' mRNA/CDS in refGene.YH.gff!!!
	# [ start, end, strand, NM, Gene, scaffold, CDS/mRNA, info ]
	push ( @{$cds{$tmp[0]}}, [ @tmp[2..6, 0, 1, 7] ]    ) if ( ($tmp[1] eq "CDS" ) or ($tmp[1] eq "mRNA") );
	$gene{$key}[0] = 0;
	$gene{$key}[1] = 0;
}
close ( I );

################################################################################
for my $key ( keys %cds ) {

	@{$cds{$key}} = sort { $a->[0] <=> $a->[0] } @{$cds{$key}};
	$index{$key}  = 0;
}
################################################################################

my ( $refStart, $refEnd, $flag, $ovlength, $totalLength );
my $lineNum = 0;
my $readNum = 0;

my $WIN_NUM    = 100;
my $totalReads = 0;
my %distribu; # 

my @bamfileT     = GetFileList( $bamfileT     );
#my @bamfileHomo  = GetFileList( $bamfileHomo  );
#my @bamfileAmb   = GetFileList( $bamfileAmb   );
#my @bamfileOther = GetFileList( $bamfileOther );

Calc ( $scaffoldId, \$totalReads, 1, @bamfileT     );
#Calc ( $scaffoldId, \$totalReads, 1, @bamfileHomo  );
#Calc ( 'all',       \$totalReads, 0, @bamfileAmb   );
#Calc ( 'all',       \$totalReads, 0, @bamfileOther );

for my $key ( keys %cvgReg ) {

	$gene{$key}[0] = MergeLen( $cvgReg{$key} );
}

print "#Id\tGene\tCDS_length\tExp_reads\tRPKM\tTotal_Exp_reads\n";
for my $key ( keys %gene ) {

	next if ( !$gene{$key}[0] ); # Means no coding region in this gene!
	print join "\t", ( (split ( /:/, $key)), @{$gene{$key}}, rpkm( @{$gene{$key}}, $totalReads ), $totalReads ) , "\n";
}
print STDERR "[Reads_number]\t$readNum\n";

for my $key ( keys %distribu ) {

	my ( $s, $g ) = split /:/, $key;
	for ( my $i = 0; $i < @{$distribu{$key}}; ++$i ) {

		print STDERR "\t", ( "#$s\t$g", $i + 1, $distribu{$key}[$i] ), "\n";
	}
}
#####################################################################################
#####################################################################################
sub MergeLen {

	my ( $region ) = @_;
	@$region = sort { $a->[0] <=> $b->[0] } @$region;

	my $preStart = $$region[0]->[0];
	my $preEnd   = $$region[0]->[1];
	my $len      = $preEnd - $preStart + 1;
	for ( my $i = 1; $i < @$region; ++$i ) {

		if ( $$region[$i]->[0] <= $preEnd && $$region[$i]->[1] > $preEnd ) {
			$len    += ( $$region[$i]->[1] - $preEnd );
			$preEnd  = $$region[$i]->[1]; 
		} elsif ( $$region[$i]->[0] > $preEnd ) {

			$len     += ( $$region[$i]->[1] - $$region[$i]->[0] + 1 );
			$preStart = $$region[$i]->[0];
			$preEnd   = $$region[$i]->[1];
		} else {

			if ( $preEnd < $$region[$i]->[1] || $preStart > $$region[$i]->[0] ) {
				die "[ERROR] Find error. You may not sorted the array by position before.\n[$preStart, $preEnd ]\n[$$region[$i]->[0], $$region[$i]->[1]]\n";
			}
		}
	}

	return $len;
}

sub GetFileList {

	my ( $fileList ) = @_;
	my @file;
	open I, $fileList or die "Cannot open file $fileList: $!\n";
	while ( <I> ) { chomp; push @file, $_; }
	close I;

	return @file;
}

sub Calc {

	my ( $scaffoldId, $totalReads, $need, @bamfile ) = @_;

	my @tmp;
	for my $bamfile ( @bamfile ) {

		open I, "/ifs1/ST_EPI/USER/huangshujia/ifs2/bin/samtools-0.1.9/samtools view -X $bamfile | " or die "Cannot open file : $bamfile\n";
		while ( <I> ) {
# FC618FNAAXX:7:94:1028:462#0     163     scaffold20      2113244 255     75M     =       2113668 499     GCCTTTT
			chomp;
			++$lineNum;
			@tmp      = split;
			next if (($tmp[1] =~ m/s/) || ($tmp[4] < 30)); # Ignore low mapQ reads and un-uniq mapping reads.
			++$readNum;
			++$$totalReads;
			next if ( ! defined $index{$tmp[2]} );

			$refStart = $tmp[3];
			$refEnd   = getEndPosition( $refStart, $tmp[5] );


			$flag = 1;
## [ start, end, strand, NM, Gene, scaffold, CDS/mRNA, info]
			my %set;
			for ( my $i = $index{$tmp[2]}; $i < @{ $cds{$tmp[2]} }; ++$i ) {
				last if $refEnd   < $cds{$tmp[2]}[$i][0];
				next if $refStart > $cds{$tmp[2]}[$i][1];

				my $key = $cds{$tmp[2]}[$i][5].":".$cds{$tmp[2]}[$i][4];
				push ( @{ $cvgReg{$key} }, [OverlapCDSregion( $refStart, $refEnd, $cds{$tmp[2]}[$i][0], $cds{$tmp[2]}[$i][1] )] );

				if ( ($scaffoldId eq $cds{$tmp[2]}[$i][5]) && !exists $set{$key} && $need ) {
					++$gene{$key}->[1]; # reads number cover this gene
				}

# reads depth distribution along the reads
				if ( $cds{$tmp[2]}[$i]->[-2] eq 'mRNA' ) {

					my @winLen = decided_win_length ( $cds{$tmp[2]}[$i][0], $cds{$tmp[2]}[$i][1]    , $WIN_NUM ); # 100 windows
					my $lowwin = decided_win ( $cds{$tmp[2]}[$i][0], $cds{$tmp[2]}[$i][1], $refStart, \@winLen ); # leftest windows
					my $larwin = decided_win ( $cds{$tmp[2]}[$i][0], $cds{$tmp[2]}[$i][1], $refEnd  , \@winLen ); # rightest windows

					if ( $cds{$tmp[2]}[$i][2] eq '-' ) { # Orientation
						$lowwin = $WIN_NUM - $lowwin - 1;
						$larwin = $WIN_NUM - $larwin - 1; 

						if ( $larwin < $lowwin ) { my $f = $lowwin; $lowwin = $larwin; $larwin = $f; }
					}

					if ( !exists $distribu{$key} ){ @{$distribu{$key}} = ( 0 ) x $WIN_NUM; }
					for ( my $i = $lowwin; $i <= $larwin; ++$i ) {
						++$distribu{$key}[$i];
					}
				}

				if ( $flag ) {
					$index{$tmp[2]} = $i if ( $flag );
					$flag           = 0;
				}
			}
			print STDERR "Have Dealed $lineNum lines.\n" if ( $lineNum % 100000 == 0 );
		}
		close ( I );
	}

	die "[ERROR] Total reads is 0!! " if ( $$totalReads == 0 );

	for my $key ( keys %cds ) {

		$index{$key}  = 0;
	}
}
#####################################################################################

sub OverlapCDSregion {
# Reture the overlap region.
    my ( $start1, $end1, $start2, $end2 ) = @_;
    my ( $start, $end ) = ( 0, 0 );

    if ( $end1 >= $start2 and $end1 <= $end2 and $start1 <= $start2 ) {
        $start = $start2;
        $end   = $end1;
    }
    elsif ( $end1 >= $start2 and $end1 <= $end2 and $start1 > $start2 ) {
        die "The end position is bigger than start position : $end1 > $start1\n" if ( $start1 > $end1 );
        $start = $start1;
        $end   = $end1;
    }
    elsif ( $start1 <= $start2 and $end1 > $end2 ) {
        die "The end position is bigger than start position : $end2 > $start2\n" if ( $start2 > $end2 );
        $start = $start2;
        $end   = $end2;
    }
    elsif ( $start1 <= $end2 and $end1 > $end2 ) {
        $start = $start1;
        $end   = $end2;
    }

    return ( $start, $end );
}


sub rpkm {

	my ( $exonlen, $read, $totalReads ) = @_;
	my ($m, $k) = ( 1000000, 1000 );
	my $rpkm = $read * ( $m * $k ) / ($totalReads * $exonlen);

	return ( $rpkm ); 
}

sub getEndPosition {

	my ( $start, $cigar ) = @_;
	my $incLength = 0;

	while ( $cigar =~ m/(\d+)(\w)/g ) {
		my( $num, $mark ) = ( $1, $2);
		$incLength += $num if ( $mark eq 'M' || $mark eq 'N' || $mark eq 'D' );
	}

	return $incLength + $start - 1;
}

sub decided_win_length {
    my ( $start, $end, $bin_num ) = @_;

    die "Region ERROR\n" if ( $end < $start );

    my $length = $end - $start;
    my $length_per_bin = int ( $length / $bin_num );
    my $remainder      = $length % $bin_num;
	my @bin_length     = ( $length_per_bin ) x $bin_num;

    if ( $length < $bin_num ) {

        @bin_length[0..($length-1)] = (1) x $length;
        return @bin_length; ## Must return the array value now! Don't contiune or it'll get wrong result!
	}

    for ( my $i = 0; $i < $remainder; ++$i ) {
       	$bin_length[$i]++;
    }
    return @bin_length;
}

sub decided_win { # Window number counts from 0;
    my ( $region_start, $region_end, $cout_site, $bin_length ) = @_;
    my $delta = $cout_site - $region_start;

    my $length = 0;
    my $i;

	#######
	for ( $i = 0; $i < @$bin_length ; ++$i ) {

		if ( $length < $delta ) {
			$length += $$bin_length[$i];
		} else {
			last;
		}
	}
	my $win = ( $i > 0 ) ? $i - 1 : 0;
	#######

	die "Out of region\nRegion: $region_start - $region_end\nCout file site: $cout_site\n" if ( !defined $win);
	return $win;
}







