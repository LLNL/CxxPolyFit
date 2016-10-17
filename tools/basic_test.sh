#! /bin/bash
rm -f output.dat
touch timings
date >> timings
echo $1 >> timings
$1 >> output.dat
diff output.dat $2/reference_output/output.dat
