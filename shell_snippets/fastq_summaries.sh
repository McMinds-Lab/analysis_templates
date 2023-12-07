## average sequence length from compressed fastq
## decompress, filter to just header and following line (getting rid of redundant qscores for this purpose), then for every other line (skipping headers) sum total length and divide by number of records.
## to get a quick estimate from top X reads, just pipe the initial zcat output through 'head -X'
zcat file.fastq.gz | grep '^@' --no-group-separator -A 1 | awk '!(NR%2) { total += length($0) } END { print(2*total/NR) }'

## longest sequence length in compressed fastq
zcat file.fastq.gz | grep '^@' --no-group-separator -A 1 | awk '!(NR%2) { currlen = length($0); if(currlen > longest) longest = currlen} END { print(longest) }'

## on Mac
zcat < file.fastq.gz | grep '^@' -A 1 | grep -v -- "^--$" | awk '!(NR%2) { currlen = length($0); if(currlen > longest) longest = currlen} END { print(longest) }'

## print seq instead of length, cat multiple files on mac
cat *.fastq.gz | zcat | grep '^@' -A 1 | grep -v -- "^--$" | awk '!(NR%2) { currlen = length($0); if(currlen > longest) {longest = currlen; seq=$0}} END { print(seq) }'


## filter fasta to only keep records whose first header field matches names in a list
## NR == FNR means execute only when the total record number is same as the per-file record number (e.g. only do these commands that store names in an array for the first file)
## RS='\n>' says to process the fasta recognizing that records are separated by a newline plus carat. Within each record, fields will be separated by whitespace including spaces and other newlines, but the first field should always be the name of the sequence
## sub and print ">" take care of formatting the newlines and carats in the output
awk 'NR == FNR {
       a[$1]++
       next
     }
     NR > FNR && FNR == 1 {
       sub(/^>/, "")
     }
     $1 in a {
       sub(/\n$/, "")
       print ">"$0
     }' seqlist.txt RS='\n>' file.fasta > file_filt.fasta
