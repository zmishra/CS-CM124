Haplophaser2 uses expectation-maximization for haplotype phasing. To run the code, first navigate to the directory with Haplophaser2 and the input file, then use the command listed below. If the output file does not exist, it will be created. The recommended partition size is 13 or 14.

Command line commands to run Haplophaser2 using Octave:
octave --no-gui --eval "Haplophaser2 (<input file>, <output file>, <partition size>)"

For example:
octave --no-gui --eval "Haplophaser2 ('toy.txt', 'toy_sol.txt', 5)"

Example input (input text file would contain the following):
2 1 2
1 0 2
1 0 2
1 0 2
0 1 1
2 1 2
1 0 2
1 0 2
1 0 2
0 1 2

Example output (output text file would contain the following):
1 1 1 0 1 1
1 0 0 0 1 1
1 0 0 0 1 1
1 0 0 0 1 1
0 0 0 1 1 0
1 1 1 0 1 1
1 0 0 0 1 1
1 0 0 0 1 1
1 0 0 0 1 1
0 0 0 1 1 1