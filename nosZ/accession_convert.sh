#download *unaligned* fungene database of seed sequences save as nosZ_unaligned_nucleotide_seqs.fa
#get list of accession numbers to be search in NCBI database, remove >
grep -ow '>\w*' nosZ_unaligned_nucleotide_seqs.fa|  sed  's/[>]//g' > accession.txt 

#pull taxonomy data from NCBI save as file with accession. <tab> taxid;King;phykum;class;order;family;genus
for acc in `cat accession.txt` ; do 
    tax=$(esearch -db nuccore -query $acc | elink -db nuccore -target taxonomy | efetch -format xml | xtract -pattern Taxon -tab ";" -first TaxId \
  -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
  -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
  -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
  -block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
  -block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
  -block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
  -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
  -group Taxon -tab ";" -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS") ;\
  echo -e "$acc\t$tax";
done > taxonomy_search
#remove taxID and add in Bacteria as Kingdom
cat taxonomy_search | sed 's/[0-9]*;//' | sed 's/'-'/'Bacteria'/'  > conversion.txt

#take fasta database and remove all nonsense except the accession number in the header
sed '/^>/s/^>\([^ ]*\) .*/>\1 /' nosZ_unaligned_nucleotide_seqs.fa > nosZ_db_accno.fa

#match header of fasta with taxonomy which outputs Kingdom;phy;etc accession number
#remove accession number and save as nosZ_db.fa

 awk 'FNR==NR{
  a[">"$1]=$2;next
}
$1 in a{
  sub(/>/,">"a[$1]" ",$1)
}1' conversion.txt nosZ_db_accno.fa |
sed '/^>/s/^>\([^ ]*\) .*/>\1 /' > nosZ_db.fa



  