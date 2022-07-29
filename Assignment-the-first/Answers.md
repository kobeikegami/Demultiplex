# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | ASCII |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | ASCII |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | ASCII |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | ASCII |
ASCII phred encoding is +33 beacuse there are some quality score characters that are small enough when converted that are too small for the +66 conversion

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. A good quality score cutoff for index/barcode reads is a quality score of 20 because the indexes are short and are required for demultiplexing our sample with high accuracy. If the indexes aren't highly accurate there will be problems with downstream analysis and sequence alignment. Using a Q score of 20 we can calculate the probability of an incorrect basecall due to chance. Then this can be mulitplied by the probabilities of the number of bases in the index called incorrectly required to reach the hamming distance of 4 for the indices (number of changes an index must undergo to be identical to another index). these independent probabilities are then multiplies to the chances of getting the other 4 indices called correctly. This end products is multiplied by the 8 choose 4 combination to yield a product that gives the overall probability of confusing an index for another index given a certain quality score. A score of 20 would result in arounf 3000 indexes of 364 million being misread. A good quality score for biological reads is less important, these reds just need to align to our reference genome and whether or not they align is more important than the mean q score at each position. Alignment is more more importance. This quality score could be much lower, lke 5
    3. Index 1, undetermined base calls:
    command: $ zcat 1294_S1_L008_R2_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
    answer: 3976613
    Index 2, undetermined base calls:
    command: $ zcat 1294_S1_L008_R3_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
    answer: 3328051

    
## Part 2
Assignment-the-first Pseudocode

1. The problem:
    Samples in a multiplexed sequencing job need to be demultiplexed (differentiated by index into respective sample types).
    Indx hopping needs to be accounted for
    Paired end reads need to be matched

2. Informative Output:
    3 file types containing numbers of the following:
        Matched paired end reads with matching indexes ( index and its rev comp)
        Indexes that are not in the barcode/index list, indexes with undetermined base calls (N's in the index base sequence)
        Unmatched index pairs

3. Examples (Unit Tests):
    Please see input FASTQ files inside the TEST-input_FASTQ directory
    Please see output files in side the TEST-output_FASTQ directory

4. Algorithm Using Pseudocode:

    -For each input file:
        -Read all files in increments of 4 lines to denote individual records
        -Go record by record and examine the sequence line of r2 and r3 inputs in order to sort them by type
            -Unknown output: record sets with poor quality indexes (contain an N for undetermined base call) have their r1 and r4 record put into two matching files of a single file type
                -r1 and r4 seq headers are modified to include the r2 index and the rev comp of the r3 index
            -Unkown output: record sets with an r2 index that is not found within the index list have their r1 and r4 record put into two matching files of a single file type (same as above)
                -r1 and r4 seq headers are modified to include the r2 index and the rev comp of the r3 index
            Hopped output: record sets with the correct r2 index but a non-matching r3 index have have their r1 and r4 record put into two matching files of a single file type (different from above)
                -r1 and r4 seq headers are modified to include the r2 index and the rev comp of the r3 index
            Matched output: record sets with an r2 index in the index list and a matching r3 have their r1 and r4 record put into two matching files of a single file type (different from above)
                -r1 and r4 seq headers are modified to include the r2 index and the rev comp of the r3 index 

5. High Level Functions:

    Reverse compliment generator
        ```takes index of the r3 record and generates its revers compliment by reversing the index character order and then generating the complimentary base sequence```
        variables: 
            index - the 2nd line of every record
            bases - "ACGTN"
        parameters:
            index length - 8
        Test: 
            provide a sample r3 file and have it output a reverse compliment of the index and compare to the index of the test file
        Return statement: should return the sequence and its reverse compliment so they can be compared side-by-side

    Phredscore converter (already built convertphred)
        ```converts ASCII score in 4th line of a record to its corresponding Q score)```
        phredscore is the 4th line in the record
        parameters:
             are range of ASCII characters and their respective q score values
        Test: 
            input: I output: 40
        Return statement: should return the base position and the q score at that position

    Output writing Functions
        ```categorizes index pairs in the 3 file types and writes them to the 2 correspionding files of each type. seq header lines modified to contain the orginal r2 record index and rev comp of r3 index
        modified r1 record put into a file and modified r4 record put into its complimentary file```
        variables:
            -r1 record
            -r4 record
            -r2 index
            -r3 index rev comp
        test:
            provide 4 unit test files with known epected output files and write output to new files and compare the file contents with CLI commands
        Return statement: should write the altered r1 and r4 sequences to paired output files of a single type


    Function to parse 4 files in parallel by record
        ```reads each of the 4 input files simultaneously in incriments of each 4 line record so that index and read pairs are group together for the downstream analyses```
        variables:
            -r1 file
            -r2 file
            -r3 file
            -r4 file
        Test:
            ???
        Return statement: should return the a single corresponding record (group of 4 lines) from each file



