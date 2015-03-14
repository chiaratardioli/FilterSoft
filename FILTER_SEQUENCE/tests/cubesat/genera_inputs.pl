#!/usr/bin/perl -w

$norb=$ARGV[0];
# We have to put the number of orbits: ./genera_inputs <norb>

open(INF,"<allorbitsSD.fla");
$i=1;
while($riga=<INF>){
    @pezzi=split(/\s+/,$riga);
    $a[$i]=$pezzi[1];
    $ecc[$i]=$pezzi[2];
    $inc[$i]=$pezzi[3];
    $Omnod[$i]=$pezzi[4];
    $omega[$i]=$pezzi[5];
    $l[$i]=$pezzi[6];
    $i++;
}
close(INF);
$max=$i-1;
#print "$max\n;

#generate a random value for the ratio A/m between 0.01 and 20
srand(time);
my $range = 19.99;
my $minimum = 0.01;

$i=1;
while($i<=$norb && $i<=$max){
#    if($i>$max)exit 1;
    open(INF,"<input.in_orig");
    open(OUTF,">input.in.$i");
    $count=0;
    while($riga=<INF>){
	$count++;
        if($count==56){
            my $random_number = rand($range) + $minimum;
            my $random_Am =sprintf("%3.3f",$random_number);
            print OUTF "Am= $random_Am\n";
        }elsif($count==72){
	    print OUTF "ic(1)= $a[$i]\n";
	}elsif($count==73){
	    print OUTF "ic(2)= $ecc[$i]\n";
	}elsif($count==74){
	    print OUTF "ic(3)= $inc[$i]\n";
	}elsif($count==76){
	    print OUTF "ic(5)= $Omnod[$i]\n";
	}elsif($count==77){
	    print OUTF "ic(6)= $omega[$i]\n";
	}elsif($count==78){
	    print OUTF "ic(7)= $l[$i]\n";
	}else{
	    print OUTF "$riga";
	}
    }
    close(OUTF);
    close(INF);
    $i++;
}
