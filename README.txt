Dev: Pietro Sette

================================================
IMPORTANT NOTES
================================================

THE PROGRAM INPUT MUST BE A DIRECTORY, AND  NOT THE DIRECT FASTA FILE. IT IS HIGHLY 
RECCOMONED TO MAKE A DIRECTORY WITH ALL OF THE SEGMENTS AND HAVE THE PATH TO THE
DIRECTORY BE THE INPUT
=====================
PROGRAM WILL WRITE FILES TITLED "ET" FOLLOWED BY THE INPUT NAME IN THE SRC DIRECTORY
Example:
Segement1.fasta(input)  ---> ETsegment1.fasta (output)

================================================
IRD DOWNLOAD INSTRUCTIONS
================================================
Select a segment
Select the option to include complete segments only
Exclude sequences after 2007
For subtype,  type in H1N1
Exclude all pH1N1 sequences
Add to a workbench

In a separate search, find the corresponding segment of A/California/04/2009 
Accession Number: FJ966082

Do a BLAST analysis with the corresponding segment of FJ966082 as the input 
Sequence and the workbench created above as the selected database. Under ouput
format, select 1000 results. Run the BLAST. When it is done select save analysis.

Go to workbench and select all of the 1000 results in the BLAST Report. Add to 
working set and then download the working set. For Download type, select Segment 
FASTA and then choose Custom Format. Select the following for header information:
Strain Name, Accession Number, Date, Geographic Location, Host

Note the location downloaded to as that will be the input of the program.

================================================
Running the Program
================================================

Command Line
=====================
Enter the top level src directory

cd src/
javac -cp . org/jcvi/psette/*.java
java -cp . org/jcvi/psette/EvolTraj

After that command the program will run. The first line should be as follows

Enter the name of the BLAST output directory:

Sample input: /Users/pietrosette/Downloads/Segments/  


Eclipse
=====================
Run with deafault configurations

================================================
Brief Program Outline
================================================


Consists of 6 classes: EvolTraj, TrajectoryCost, ETstrain, DistanceComparator,
YearComparator, and MacroData
Order:
EvolTraj ->  TrajectoryCost -> ETstrain -> TrajectoryCost -> DistanceComparator ->
TrajectoryCost -> YearComparator -> TrajectoryCost -> Macrodata -> EvolTraj

EvolTraj
=====================
Methods included: readSequenceFromFile, getHeader, getSequence, size (All self explanatory)
Functions included: main, arrayfiller,psalgo, psalgo2, print

Main reads in the directory, finds the fasta files, and calls the other functions

Arrayfiller creates and fills the arrays of headers and begining/end positions

psAlgo reads in information from the sequence and the header and creates values for
genetic difference (from the refrence strain) and difference in year (from the 
refrence strain). From this,the algorithm finds the lowest slope from the 
reference point(which is centered around the origin). After it finds the point 
with the lowest slope, it re-centers about that point. The algorithm re-does the 
calculation so that the x and y coordinates are years from the point and genetic
differences from the point. The point is at the origin. It does this as needed 
until the rightmost point is reached. It creates a lowerbounded convex hull

psAlgo2 uses this estimated mutation rate and creates a region bounded by this 
mutation rate and the convex hull. It adds all of the x&y coordinates of the points
in that region and estimates a new, more precise and more accurate mutation rate

Print is the print function just prints out the psAlgo2 mutation rate and the 
summary (elaborated on in TrajectoryCost).


TrajectoryCost
=====================
Functions included: initilaizer, comparing, interp, and finaldata

initalizer initializes the datastructure array of strains of type ETstrain

comparing does the penalty computation by inserting each sequence one at a time 
and creating a multiple sequence alignment, finding the point where adding more
sequences is no longer beneficial. It adds the sequence by distance from the 
psAlso2 estimated mutation rate (IMPORTANT: TrajectoryCost.intializer(ads, 
macroarray, header, sequence, startPos, endPos, smg, bigcounter, u, psa_results,
high, medium); will sort by distance from psAlgo2. TrajectoryCost.intializer(
, macroarray, header, sequence, startPos, endPos, fmg, bigcounter, u, 
psa_results, high, medium); will sort by distance from lowest point. This line 
of code is located at line 132 in EvolTraj, in the main function). Comparing 
also creates a table of each mutation, storing the position, the year, from what 
nucleotide and to what nucleotide.

interp interprets the data, finding  the most common geographic location and host 
in the list of ancestral candidates. This also prints out the aforementioned table

finaldata finds the most common geographic location and host for all segments, 
determining an evolved from location and host (summary).


ETstrain
=====================
Includes the variables: distancefrom, ycord, seq, strain_name, accesion_number, 
geographic_location, host, year, xcord

DistanceComparator
=====================
Sorts the sequences based on distance from the input line

YearComparator
=====================
Sorts the sequence based on year

MacroData
=====================
Assembles all the data used for finaldata
