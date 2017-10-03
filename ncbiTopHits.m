function [accList] = ncbiTopHits (accInput, N)
% [AccessionList] = ncbiTopHits[AccessionInput, N]
% AccessionInput defines the initial Accession number of the sequence you
% wish to blast and N defines the number of hits that it will return. If N
% is greater than the total number of Hits, on the top hits will be
% returned. 
disp('Initializing Top NCBI Hits...')
gb_data = getgenbank(accInput);
seqlng = length(gb_data.Sequence);
if seqlng >= 200
    seqlng = 200;
end

blstseq = gb_data.Sequence(1:seqlng);
[reqID,reqTime] = blastncbi(blstseq,'blastn');
fprintf('Attempting to get Blast Data from NCBI. \nEstimated time: %f sec \n',reqTime)
blast_data = getblast(reqID,'WaitTime',reqTime);

if  N > length(blast_data.Hits)
    N = length(blast_data.Hits);
    fprintf('Requested number of hits is too large. \n')
    fprintf('List will now contain only  %g hits. \n',N)
end

for ii = 1:N
    segarray = strfind(blast_data.Hits(ii).Name,'|');
    j = segarray(1,3);
    k = segarray(1,4);
    x = blast_data.Hits(ii).Name(1,j+1:k-3);
    accList{ii,:} = x;
end

end
% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits.