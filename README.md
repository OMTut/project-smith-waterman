# Our Project Repo!

s3_it is a fully working iterative implementation of the Smith-Waterman algorithm used to find local alignments.

Because the size of the genomes are huge, I'm only printing out an alignment scores. The apps all take fasta files from the /dataset directory.

### Compiling
```
g++ -o whatever s3_it.cpp
```
### Running
```
./whatever ../dataset/filename.fasta ../dataset/filename2.fasta
```

### Fasta Files
They can be edited just like txt files. Just open it up, make your changes, and save.

## Smith-Waterman Algo
A high level (because that's all I understand at this point) is that the SW algo creates a matrix from the two sequences and fills it with scores. The scores represent how close the elements in the sequences match. Each score is computed using a "match score", "mismatch_penalty", and "gap_penalty".

Scoring<BR>
In our case, we're scoring match +=2, mismatch_penalty -=1, and gap_penalty -= 1. A match ends up being 2pts and everything else is 0.

Example: seq1='GGG' and seq2='GGG'<BR>
Alginment Score: 6

Example: seq1='GGG' and seq2='TTT'<BR>
Alignment Score: 0

Example seq1='GGG" and seq2='TTA'<BR>
ALignment Score: 4

Note: in the full genome, we are diving the millions of characters into blocks of 10000. I'm not 100% sure, but a score of 7500 out of 10000 means that it's relatively similar. I'm not familiar enough with bioinformatics to know exactly how or why yet.