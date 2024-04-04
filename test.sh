#!/bin/bash

total_sequential_time=0
total_unrolling_time=0

program_name=gepp_omp

for i in {1..10}
do
    # output=$(./$program_name 1000)
    output=$(./$program_name 1000 4)
    
    sequential_time_value=$(echo "$output" | grep "sequential calculation time" | awk '{print $4}')

	if ! [[ $sequential_time_value =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        echo "Error: Sequential time value '$sequential_time_value' is not a number."
        exit 1
    fi
    total_sequential_time=$(echo "$total_sequential_time + $sequential_time_value" | bc)
    
    unrolling_time_value=$(echo "$output" | grep "unrolling and blocking calculation time" | awk '{print $6}')

    if ! [[ $unrolling_time_value =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        echo "Error: Unrolling time value '$unrolling_time_value' is not a number."
        exit 1
    fi
    total_unrolling_time=$(echo "$total_unrolling_time + $unrolling_time_value" | bc)
done

average_sequential_time=$(echo "$total_sequential_time / 10" | bc -l)
average_unrolling_time=$(echo "$total_unrolling_time / 10" | bc -l)

echo "Average Sequential Calculation Time: $average_sequential_time seconds"
echo "Average Unrolling and Blocking Calculation Time: $average_unrolling_time seconds"
