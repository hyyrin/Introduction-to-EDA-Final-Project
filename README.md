# README

---

## 1. Project Overview

This project explores solutions for the 2024/2025 ICCAD CAD Contest, Problem B. The goal is to optimize digital circuits for power, timing, and area by intelligently applying multibit flip-flop (MBFF) banking and debanking. We implemented and compared three distinct clustering algorithms for this task: **Greedy**, **K-Means**, and **Mean-Shift**. All algorithms are built upon a robust two-phase framework that separates logical transformation from physical placement to ensure correct, overlap-free results.

## 2. File and Directory Structure

The submission zip file contains the following structure:

```
.
├── README.md                # This file
├── report.pdf               # Final project report
├── slide.pdf                # Final presentation slides
├── greedy/
│   ├── greedy.cpp           # Source code for the Greedy algorithm
│   └── Makefile             # Makefile to compile the Greedy version
├── k-mean/
│   ├── k-mean.cpp           # Source code for the K-Means algorithm
│   └── Makefile             # Makefile to compile the K-Means version
├── mean-shift/
│   ├── MeanShift.cpp        # Source code for the Mean-Shift algorithm
│   └── Makefile             # Makefile to compile the Mean-Shift version
├── testcases/               # Directory containing all test cases
├── docs/                    # Files I reference.
└── results/                 # Directory containing detailed output logs
```

## 3. Reproduce the Results

We use the testcases provided by the ICCAD 2024 CAD Contest, Problem B.

### 3.1 Compilation

Each algorithm is located in its own directory and can be compiled independently using the provided `Makefile`.

To compile a specific algorithm (e.g., K-Means), navigate to its directory and run `make`.

**Example:**
```bash
# Navigate into the k-mean directory
cd k-mean/

# Run make to compile
make
```
This will generate an executable named `optimizer` in the current directory. The process is the same for the `greedy` and `mean-shift` directories.

### 3.2. Running a Single Algorithm

The program takes the paths to the input files and a final base name for the output files.

**Syntax:**
```bash
./optimizer inputs/<inputFile1.txt> outputs/<outputName>
```

**Example:**
To run the executable on `testcase1_0614`:
```bash
./optimizer inputs/testcase1_0614.txt outputs/testcase1_0614_output.txt
```
This command will generate the output file `testcase1_0614_output.txt` in the `outputs` directory.


### 3.3 Using Helper Scripts (Optional)
For convenience, a script `run.sh` is provided to automate the execution of an algorithm across all test cases. You could just run
```bash
./run.sh
```
All the output files would be generated in the `outputs` directory at once.

## 4. Evaluation and Verification

We use the executables provided by the ICCAD 2024 CAD Contest, Problem B for scoring and sanity checks.

### 4.1. Scoring with the Official Evaluator
To calculate the score for a generated output, use the provided evaluator named `main`.

**Syntax:**
```bash
./main inputs/<inputFile1.txt> outputs/<outputName>
```

**Example:**
To run `testcase1_0614`,
```bash
./main inputs/testcase1_0614.txt outputs/testcase1_0614_output.txt
```
This will print the performance metrics (TNS, Power, Area) and the final score.

### 4.2 Using Helper Scripts (Optional)
For convenience, a script `run_v2.sh` is provided to automate the execution of the evaluator accross all testcases. You could just run
```bash
./run_v2.sh
```
All the results would be recored in the `output.log`.


### 4.3 Verifying Pin Connections
To verify the correctness of the generated pin mappings, use the provided sanity check tool named `sanity_20240801`
**Syntax:**
```bash
./sanity_20240801 inputs/<inputFile1.txt> outputs/<outputName>
```

**Example:**
To run `testcase1_0614`,

```bash
./sanity_20240801 inputs/testcase1_0614.txt outputs/testcase1_0614_output.txt
```

---
