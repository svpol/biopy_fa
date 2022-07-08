The package consists of:
1. mat_fa - allows to perform matrix processes with nucleotide sequences. This file also contains an example usage pipeline.
2. search_fa - allows to search records having a particular string in their description.
3. sort_fa - allows to sort files by record ID or by sequence length. By now - for ungzipped files only.

The functions are designed to work with fasta or fastq files. Other formats applicable to Biopython are not forbidden explicitly, but once the user chooses any different format, the functions' behaviour may be unpredicted.

There are several files to test the functionality in test_files directory.