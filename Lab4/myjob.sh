#!/bin/bash
# Simple executable to run on the SMUHPC nodes
# Dan Reynolds
# May 2013

# start with "2"
echo "2 is prime"
Nprime=1

# loop until we've found the requisite number of primes
# Note: the command-line argument is the number primes to compute
n=2
while [ $Nprime -lt $1 ]; do

   n=$(($n + 1))
   is_prime=1

   # loop over all values smaller than n
   for ((i=2; i<$n; i++)); do
   
      # stop once we've exceeded sqrt(n)
      i2=$(($i * $i))
      if [ $i2 -gt $n ]; then
         break
      fi

      # check if n is evenly divisible by i
      md=$(( $n % $i ))
      if [ $md -eq 0 ]; then
         is_prime=0
	 break
      fi
   done

   # check for prime
   if [ $is_prime -eq 1 ]; then
       echo "$n is prime"
       Nprime=$(( $Nprime + 1 ))
   fi

done
