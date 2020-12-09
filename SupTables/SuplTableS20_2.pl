

chdir('C:\Users\CAntonescu\Desktop\try');

use warnings;

open IN, "intervals_sorted.txt";
open OUT, ">intervals_disjoint.txt";
open OUT2, ">intervals_contained.txt";
open OUT3, ">intervals_intersect.txt";

print OUT "chr\tstart\tend\n";
print OUT2 "chr\tstart1\tend1\tstart2\tend2\n";
print OUT3 "chr\tstart1\tend1\tstart2\tend2\n";

$lastchr = "";
$laststart = 0;
$lastend = 0;
while (<IN>) {
    next if $.==1;
    chomp;
    ($chr, $start, $end) = split /\t/;
    if ($chr eq $lastchr && $lastend >= $start) {
        if ($lastend >= $end) {
            print OUT2 "$chr\t$laststart\t$lastend\t$start\t$end\n";
            next;
        }
        else {
            print OUT3 "$chr\t$laststart\t$lastend\t$start\t$end\n";
            next;
        }
    }
    print OUT "$chr\t$start\t$end\n";
    $lastchr = $chr;
    $laststart = $start;
    $lastend = $end;
}
