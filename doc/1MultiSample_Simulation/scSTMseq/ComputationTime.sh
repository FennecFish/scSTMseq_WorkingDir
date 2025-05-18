#!/bin/bash
#SBATCH --job-name=MeanTimeCalculation_linearMixed
#SBATCH --output=MeanTimeCalculation_linearMixed%j.out
#SBATCH --error=MeanTimeCalculation_linearMixed%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

# Job ID
ARRAY_JOB_ID=52350824

# output file
OUTPUT_FILE="mean_time_${ARRAY_JOB_ID}.txt"
echo -e "Array Job ID\tComputing Time" > "$OUTPUT_FILE"

# calculate the total time and count of completed tasks
total_time=0
count=0

# elapsed times for COMPLETED tasks in the array job
for time in $(sacct -j "$ARRAY_JOB_ID" --state=COMPLETED -o Elapsed -n | awk '{print $1}'); do
    if [[ "$time" == *-* ]]; then
        IFS='-' read -r days time_only <<< "$time"
    else
        days=0
        time_only=$time
    fi

    # Split HH:MM:SS
    IFS=: read -r hours minutes seconds <<< "$time_only"
    
    # Strip leading zeros and handle time conversion
    days=$((10#${days#0}))
    hours=$((10#${hours#0}))
    minutes=$((10#${minutes#0}))
    seconds=$((10#${seconds#0}))
    
    # Calculate elapsed seconds including days
    elapsed_seconds=$((days * 86400 + hours * 3600 + minutes * 60 + seconds))
    
    # Add to total time and increment count
    total_time=$((total_time + elapsed_seconds))
    count=$((count + 1))
done

# Calculate the mean time if there are any completed tasks
if (( count > 0 )); then
    mean_seconds=$((total_time / count))
    
    # Calculate days, hours, minutes, and seconds from mean_seconds
    mean_days=$((mean_seconds / 86400))
    remaining_seconds=$((mean_seconds % 86400))
    mean_hours=$((remaining_seconds / 3600))
    remaining_seconds=$((remaining_seconds % 3600))
    mean_minutes=$((remaining_seconds / 60))
    mean_seconds=$((remaining_seconds % 60))
    
    # Format the mean time as D-HH:MM:SS
    mean_time_formatted=$(printf "%d-%02d:%02d:%02d" "$mean_days" "$mean_hours" "$mean_minutes" "$mean_seconds")
    echo -e "${ARRAY_JOB_ID}\t${mean_time_formatted}" >> "$OUTPUT_FILE"
    echo "Mean time for COMPLETED tasks written to $OUTPUT_FILE"
else
    echo -e "${ARRAY_JOB_ID}\tNo completed tasks found" >> "$OUTPUT_FILE"
    echo "No completed tasks found for array job $ARRAY_JOB_ID."
fi
