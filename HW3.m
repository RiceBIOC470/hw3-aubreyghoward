% GB comments
1.	100
2a. 70 Calculated incorrectly. You are using the total length of the aligned sequence which is not the same as the length of the coding sequence. You need to divide the shared number of homology by the total length of erk1 or erk2. 
2b. 70 same issue as 2a
2c. 70 same issue as 2a. 
3a 100 
3b. 100
3c. 100  	
Overall: 87


%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

%Adam Howard: Please see image saved in repository


%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

erk1 = getgenbank('NM_002746');%gets genbank data for ERK1
erk2 = getgenbank('NM_002745');%gets genbank data for ERK2

[score,align] = swalign(erk1,erk2,'Alphabet','nt','Showscore',false);
%Performs the Smith-Watterman alignment for ERK1 & ERK2. Returns the score
%and, most importantly, the alignment. Showscore is set to false for speed
%of processing.

fracmatch = count(align(2,:),'|');%counts the # of matches
fracmatch = fracmatch/length(align(3,:));%gets the fraction of matches
fprintf('NT Faction Match = %f \n',fracmatch)%prints result


% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

[score,align] = swalign(erk1.CDS.translation,erk2.CDS.translation,'Alphabet','AA','Showscore',false);
%Performs the Smith-Watterman alignment for ERK1 & ERK2. Returns the score
%and, most importantly, the alignment. Showscore is set to false for speed
%of processing.

fracmatch = count(align(2,:),'|')+count(align(2,:),':');%counts the # of matches
fracmatch = fracmatch/length(align(3,:));%gets the fraction of matches
fprintf('AA Faction Match = %f \n',fracmatch)%prints result


% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

Merk1 = getgenbank('X64605');%gets genbank data for mouse ERK1 
Merk2 = getgenbank('AK087925');%gets genbank data for mouse ERK2 
showgraph = false;%turns on graphs, if desired

[score,align] = swalign(Merk1,erk1,'Alphabet','nt','Showscore',showgraph);
fprintf('ERK1 NT Seq Comparison = %f \n',score)
for ii = 1:3
fprintf('Alignment = %s \n',align(ii,:))
end%prints the data for ERK1 NT compare
[score,align] = swalign(Merk1.CDS.text(17,15:53),erk1.CDS.translation,'Alphabet','AA','Showscore',showgraph);
fprintf('ERK1 AA Seq Comparison = %f \n',score)
for ii = 1:3
fprintf('Alignment = %s \n',align(ii,:))
end%prints the data for ERK1 AA compare 
[score,align] = swalign(Merk2,erk2,'Alphabet','nt','Showscore',showgraph);
fprintf('ERK2 NT Seq Comparison = %f \n',score)
for ii = 1:3
fprintf('Alignment = %s \n',align(ii,:))
end%prints the data for ERK2 NT compare
[score,align] = swalign(Merk2.CDS.translation,erk2.CDS.translation,'Alphabet','AA','Showscore',showgraph);
fprintf('ERK2 AA Seq Comparison = %f \n',score)
for ii = 1:3
fprintf('Alignment = %s \n',align(ii,:))
end%prints the data for ERK2 AA compare
fprintf('Part 2 Completed \n')

%%Adam Howard: First, by comparing the mouse/human match scores between the ERK1 and
%%ERK2, its clear that both ERK1 and ERK2 have strong sequence homology between two speices. 
%ERK 2 has a stronger match score in both comparisons and matches over a much
%%longer number of nuclaic acids and amino acid residues. 

%% Problem 3: using blast tools programatically
clear all
% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 
%[accList] = ncbiTopHits('NM_002745',10)

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

% [hiHum,hiXen] = bestgenematches ('NC_000017')

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

%Please note, that running this code, can take a significant amount of
%time.

[hiHum,hiXen] = bestgenematches ('U60259')%Human input gene for MUC1 transmembrane protein

%Output: hiHum ='EF670712' is mRNA from a human isoform of MUC1; hiXen =
%'XM_004026892' is the mRNA from a lowland Gorilla for MUC1. 
%%
[hiHum2,hiXen2] = bestgenematches ('X51749')%Non-human input

%Output: no human homolouges found; hiXen2 = 'X51749' is the same gene that
%was used to find in Drosophila. There are several varients for white eye 
%pigment, which is why I expect I got this resutl. 

