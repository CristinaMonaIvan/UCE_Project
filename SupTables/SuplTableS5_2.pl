



use warnings;

@f = glob "*.vcf";

open OUT, ">mut_all.txt";

print OUT "CHROM\tPOS\tEND\n";

for $f (@f) {
    open IN, $f;
    while (<IN>) {
        next if /^#/;
        ($chrom, $pos, undef, $ref) = split /\t/;
        $end = $pos+length($ref)-1;
        print OUT "$chrom\t$pos\t$end\n";
    }
}



