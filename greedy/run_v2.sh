#!/bin/bash

# clean log file
> output.log

KEYWORDS="area change to|timing changed to|power changed to|tns:|FF area|Final score:|Violated bin"

INPUT_FILES=(
  inputs/testcase1_0614.txt
  inputs/testcase1_0812.txt
  inputs/testcase2_0812.txt
  inputs/testcase3.txt
  inputs/hiddencase01.txt
  inputs/hiddencase02.txt
  inputs/hiddencase03.txt
  inputs/hiddencase04.txt
)

for input in "${INPUT_FILES[@]}"; do
    filename=$(basename "$input")
    output="outputs/${filename%.txt}_output.txt"

    echo "Processing $input ..." | tee -a output.log

    # write to log file
    ./main "$input" "$output" 2>&1 | grep -E "$KEYWORDS" >> output.log

    echo "" >> output.log
done
