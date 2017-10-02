function [hiHum,hiXen] = bestgenematches(accInput)

[accList] =  ncbiTopHits(accInput,50);
ref_data = getgenbank(accInput);
for ii = 1:length(accList)
    checkgene = getgenbank(accList{ii,1});
    if strfind(checkgene.SourceOrganism(1,:),'Homo sapiens') > 0        
        humList{ii,:} = accList{ii,1};
    else
        xenList{ii,:} = accList{ii,1};
    end
end
fprintf('Completed Blast Data Collection from NCBI.')

showgraph = false;
hiHumScore = -1;

disp('Now Processing Human gene data. Please wait...')
score = 0;  %resets the score value between human & non-human genes
if length(humList)<= 0
    disp ('No human homologues found')
else
    humList = humList(~cellfun('isempty',humList)); %Removes empty cells. Referenced to :"https://www.mathworks.com/matlabcentral/answers/27042-how-to-remove-empty-cell-array-contents"
    for ii = 1:length(humList)
    humcheck = getgenbank(humList{ii,1});
    [score] = swalign(ref_data.Sequence,humcheck.Sequence,'Alphabet','nt','Showscore',showgraph);
        if score > hiHumScore
            hiHumScore = score;
            hiHum = humList{ii,:};
        end
    disp('Still processing Human gene list. Please wait...')
    end
end
     
hiXenScore = -1;
disp('Now Processing Non-Human gene data. Please wait...')
if length(xenList)<= 0
    disp ('No non-human homologues found')
else
    xenList = xenList(~cellfun('isempty',xenList)) %Removes empty cells. See earlier reference.
    for ii = 1:length(xenList)
    xencheck = getgenbank(xenList{ii,1});
    [score] = swalign(ref_data.Sequence,xencheck.Sequence,'Alphabet','nt','Showscore',showgraph);
        if score > hiXenScore
            hiXenScore = score;
            hiXen = xenList{ii,:};
        end
    disp('Still processing Non-Human gene list. Please wait...')
    end
end
disp('bestgenematches analysis complete')
end
% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found.