#!/usr/bin/perl
use strict;
use warnings;

my %pop_samples = (
    'AC' => ['Biyu', 'H0809', 'Hongyang', 'Hort16A', 'Huangyang', 'Jinmi', 'Jinpai', 'Zps'],
    'Ae' => ['Blank', 'CY38', 'DG01', 'GZ26', 'HN01', 'LD19', 'MQ32', 'QS41', 'RY10', 'YL01', 'YX01', 'YX02'],
    'F1' => ['HH01']
);

my %pop_expected_alleles = (
    'AC' => 8 * 2,
    'Ae' => 12 * 2,
    'F1' => 1 * 2
);


my $SV_THRESHOLD = 50;

my %sample_to_pop;
foreach my $pop (keys %pop_samples) {
    foreach my $sample (@{$pop_samples{$pop}}) {
        $sample_to_pop{$sample} = $pop;
    }
}

if (@ARGV != 1) {
    print "perl vcf_freq.pl input.vcf(.gz)\n";
    print "CHROM\tPOS\tID\tREF\tALT\tAC_Freq\tAe_Freq\tF1_Freq\tAC_Count\tAe_Count\tF1_Count\tTotal_Count\tAC_MissingRate\tAe_MissingRate\tF1_MissingRate\tTotal_MissingRate\tVariant_Type\tMax_Length\tLength_Change\n";
    print "SVtype\n";
    print "SNP\n";
    print "Indel\n";
    print "SV\n";
    print "Max_Length\n";
    print "Length_Change\n";
    exit(1);
}

my $vcf_file = $ARGV[0];


my $input_handle;
if ($vcf_file =~ /\.gz$/) {
    open($input_handle, "gunzip -c $vcf_file |") or die "cannot open $vcf_file: $!";
} else {
    open($input_handle, "<", $vcf_file) or die "cannot open $vcf_file: $!";
}

print "CHROM\tPOS\tID\tREF\tALT\tAC_Freq\tAe_Freq\tF1_Freq\tAC_Count\tAe_Count\tF1_Count\tTotal_Count\tAC_MissingRate\tAe_MissingRate\tF1_MissingRate\tTotal_MissingRate\tVariant_Type\tMax_Length\tLength_Change\n";

my @header;
my %sample_index;
my $data_started = 0;

while (my $line = <$input_handle>) {
    chomp $line;
    
    if ($line =~ /^##/) {
        next;
    }
    
    if ($line =~ /^#CHROM/) {
        @header = split(/\t/, $line);
        $header[0] =~ s/^#//;
        
        for (my $i = 9; $i < @header; $i++) {
            my $sample_name = $header[$i];
            if (exists $sample_to_pop{$sample_name}) {
                $sample_index{$sample_name} = $i;
            }
        }
        $data_started = 1;
        next;
    }
    
    next unless $data_started;
    
    my @fields = split(/\t/, $line);
    my ($chrom, $pos, $id, $ref, $alt) = @fields[0..4];
    $id = '.' if $id eq '';
    

    my ($variant_type, $max_length, $length_change) = classify_variant_by_length($ref, $alt);
    

    my %pop_stats = (
        'AC' => {total => 0, alt => 0, observed => 0, missing => 0},
        'Ae' => {total => 0, alt => 0, observed => 0, missing => 0},
        'F1' => {total => 0, alt => 0, observed => 0, missing => 0}
    );
    

    foreach my $sample_name (keys %sample_index) {
        my $pop = $sample_to_pop{$sample_name};
        my $genotype_field = $fields[$sample_index{$sample_name}];
        
        if ($genotype_field eq './.' || $genotype_field eq '.') {
            $pop_stats{$pop}{missing} += 2;
            $pop_stats{$pop}{observed} += 2;
            next;
        }
        
        my @parts = split(/:/, $genotype_field);
        my $gt = $parts[0];
        my @alleles = split(/[\/\|]/, $gt);
        
        foreach my $allele (@alleles) {
            $pop_stats{$pop}{observed}++;
            
            if ($allele eq '.') {
                $pop_stats{$pop}{missing}++;
            } else {
                $pop_stats{$pop}{total}++;
                if ($allele ne '0') {
                    $pop_stats{$pop}{alt}++;
                }
            }
        }
    }
    
    my @freqs;
    my @counts;
    my @missing_rates;
    my $total_alt = 0;
    my $total_alleles = 0;
    my $total_observed = 0;
    my $total_expected = 0;
    my $total_missing = 0;
    
    foreach my $pop (qw(AC Ae F1)) {
        my $stats = $pop_stats{$pop};
        my $expected = $pop_expected_alleles{$pop};
        
        my $freq = $stats->{total} > 0 ? 
                   sprintf("%.6f", $stats->{alt} / $stats->{total}) : "0.000000";
        
        my $count_str = "$stats->{alt}/$stats->{total}";
        
        my $missing_rate = "0.0000";
        if ($expected > 0) {
            $missing_rate = sprintf("%.4f", $stats->{missing} / $expected);
        }
        
        push @freqs, $freq;
        push @counts, $count_str;
        push @missing_rates, $missing_rate;
        
        $total_alt += $stats->{alt};
        $total_alleles += $stats->{total};
        $total_observed += $stats->{observed};
        $total_expected += $expected;
        $total_missing += $stats->{missing};
    }
    
    my $total_missing_rate = "0.0000";
    if ($total_expected > 0) {
        $total_missing_rate = sprintf("%.4f", $total_missing / $total_expected);
    }
    
    print join("\t", 
        $chrom, $pos, $id, $ref, $alt,
        @freqs,
        @counts,
        "$total_alt/$total_alleles",
        @missing_rates,
        $total_missing_rate,
        $variant_type,
        $max_length,
        $length_change
    ), "\n";
}

close($input_handle);

unless ($data_started) {
    print STDERR "can not find var\n";
}


sub classify_variant_by_length {
    my ($ref, $alt_str) = @_;
    
    my @alts = split(/,/, $alt_str);
    my $ref_len = length($ref);
    

    my $max_length = $ref_len;
    my $max_length_change = 0;
    

    my @lengths = ($ref_len);
    

    foreach my $alt (@alts) {
        my $alt_len = length($alt);
        push @lengths, $alt_len;
        
        $max_length = $alt_len if $alt_len > $max_length;
        
        my $len_diff = abs($alt_len - $ref_len);
        $max_length_change = $len_diff if $len_diff > $max_length_change;
    }
    
    for (my $i = 0; $i < @alts; $i++) {
        for (my $j = $i + 1; $j < @alts; $j++) {
            my $len_diff = abs(length($alts[$i]) - length($alts[$j]));
            $max_length_change = $len_diff if $len_diff > $max_length_change;
        }
    }
    
    my $all_length_one = 1;
    foreach my $len (@lengths) {
        if ($len != 1) {
            $all_length_one = 0;
            last;
        }
    }
    
    if ($all_length_one) {
        return ("SNP", 1, 0);
    }
    
    if ($max_length >= $SV_THRESHOLD) {
        return ("SV", $max_length, $max_length_change);
    } else {
        return ("Indel", $max_length, $max_length_change);
    }
}