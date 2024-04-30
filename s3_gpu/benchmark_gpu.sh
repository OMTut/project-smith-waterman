#!/bin/bash

# Name of the C++ executable
EXECUTABLE="./app"
ARGS="../dataset/NC_000913_3_trimmed_1_6.fasta ../dataset/NC_002695_2_trimmed_3k.fasta"

# Output file
OUTPUT_FILE="benchmark_results_gpu.csv"
echo "Run No, Real Time, User Time, Sys Time" > $OUTPUT_FILE

# Clear the output file
> $OUTPUT_FILE

# Array to store times
num_runs=100
total_real_time=0
total_user_time=0
total_sys_time=0

# Set the TIMEFORMAT variable
TIMEFORMAT='%3R %3U %3S'

# Run the the program and measure 10 times
for (( i=0; i<100; i++ )); do
   echo "Currently on run $i"
   # Capture the real, user, and syst time for each run
   { time $EXECUTABLE $ARGS > program_output.txt; } 2> time_output.txt
   read -r real_time user_time sys_time < time_output.txt
   output=$(cat program_output.txt)

   # Extract the alignment score from the output
   alignment_score=$(echo "$output" | grep -oP 'Alignment score: \K\d+')
   
   # Echo out the individual and total results
   echo "Run $i : real $real_time, user $user_time, sys $sys_time"
   total_real_time=$(echo "$total_real_time + $real_time" | bc)
   total_user_time=$(echo "$total_user_time + $user_time" | bc)
   total_sys_time=$(echo "$total_sys_time + $sys_time" | bc)

   # Writ the output to a file
   echo "$i, $real_time, $user_time, $sys_time" >> $OUTPUT_FILE
done
average_real_time=$(echo "scale=2; $total_real_time / $num_runs" | bc)
average_user_time=$(echo "scale=2; $total_user_time / $num_runs" | bc)
average_sys_time=$(echo "scale=2; $total_sys_time / $num_runs" | bc)
echo "Average Times:"
echo "real time: $average_real_time, user time: $average_user_time, sys time: $average_sys_time"

