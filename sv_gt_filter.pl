#!/usr/bin/perl
use strict;
use warnings;

my %pop1_samples = (
    'Biyu' => 1,
    'H0809' => 1,
    'Hongyang' => 1,
    'Hort16A' => 1,
    'Huangyang' => 1,
    'Jinmi' => 1,
    'Jinpai' => 1,
    'Zps' => 1
);

my %hybrid_samples = (
    'HH01' => 1
);

my $vcf_file = shift @ARGV or die$!;

open my $in, '<', $vcf_file or die "cannot: $vcf_file\n";

my @original_header_samples;
my %sample_indices;
my @reordered_samples;
my @new_header_fields;

while (<$in>) {
    chomp;
    
    if (/^##/) {
        print "$_\n";
        next;
    }
    
    if (/^#CHROM/) {
        my @fields = split /\t/;
        
        @new_header_fields = @fields[0..8];
        
        for my $i (9..$#fields) {
            $sample_indices{$fields[$i]} = $i;
            push @original_header_samples, $fields[$i];
        }
        

        @reordered_samples = ();
        

        foreach my $sample (@original_header_samples) {
            if (exists $pop1_samples{$sample}) {
                push @reordered_samples, $sample;
            }
        }

        foreach my $sample (@original_header_samples) {
            if (!exists $pop1_samples{$sample} && !exists $hybrid_samples{$sample}) {
                push @reordered_samples, $sample;
            }
        }
        
        foreach my $sample (@original_header_samples) {
            if (exists $hybrid_samples{$sample}) {
                push @reordered_samples, $sample;
            }
        }
        
        print join("\t", @new_header_fields, @reordered_samples), "\n";
        next;
    }
    
    my @fields = split /\t/;
    
    my @fixed_fields = @fields[0..8];
    
    my %pop1_haplotypes; 
    my %pop2_haplotypes; 
    
    for my $sample (@original_header_samples) {
        my $idx = $sample_indices{$sample};
        my $gt_field = $fields[$idx];
        
        if ($gt_field ne '.|.' && $gt_field ne './.') {
            my ($hap1, $hap2) = split /[|\/]/, $gt_field;
            
            if (exists $pop1_samples{$sample}) {
                if ($hap1 ne '.') { $pop1_haplotypes{$hap1}++; }
                if ($hap2 && $hap2 ne '.') { $pop1_haplotypes{$hap2}++; }
            } elsif (!exists $hybrid_samples{$sample}) {
                if ($hap1 ne '.') { $pop2_haplotypes{$hap1}++; }
                if ($hap2 && $hap2 ne '.') { $pop2_haplotypes{$hap2}++; }
            }
        }
    }
    
    next if scalar(keys %pop1_haplotypes) == 0 || scalar(keys %pop2_haplotypes) == 0;
    
    my $completely_different = 1;
    
    foreach my $hap1 (keys %pop1_haplotypes) {
        if (exists $pop2_haplotypes{$hap1}) {
            $completely_different = 0;
            last;
        }
    }
    

    if ($completely_different) {
        my @reordered_genotypes;
        foreach my $sample (@reordered_samples) {
            my $idx = $sample_indices{$sample};
            push @reordered_genotypes, $fields[$idx];
        }
        print join("\t", @fixed_fields, @reordered_genotypes), "\n";
        

    }
}

close $in;



my @maohua_samples;
foreach my $sample (@original_header_samples) {
    if (!exists $pop1_samples{$sample} && !exists $hybrid_samples{$sample}) {
        push @maohua_samples, $sample;
    }
}

for (my $i = 0; $i < @reordered_samples; $i++) {
    print STDERR "  " . ($i+1) . ". $reordered_samples[$i]\n";
}