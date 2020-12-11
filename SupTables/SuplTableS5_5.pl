
use warnings;

open IN, "morethanone.txt";

while (<IN>) {
    next if $.==1;
    chomp;
    ($chrom, $pos) = split /\t/;
    $dupl{$chrom}->{$pos} = 1;
}

@f = glob "*.vcf";

open OUT, ">mutevents.txt";
$header = 1;

for $f (@f) {
    open IN, $f;
    while (<IN>) {
        next if /^##/;
        if (/^#/) {
            if ($header) {
                s/^#//;
                print OUT;
                $header = 0;
            }
            next;
        }
        ($chrom, $pos) = split /\t/;
        print OUT unless $dupl{$chrom}->{$pos};
    }
}

