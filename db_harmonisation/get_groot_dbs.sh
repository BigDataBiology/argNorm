mkdir resfinder && cd $_
wget https://bitbucket.org/genomicepidemiology/resfinder_db/get/dc33e2f9ec2c.zip -O resfinder.zip
unzip resfinder.zip
awk 'FNR==1{print ""}1' genomic*/*.fsa > resfinder-refs.fna
awk '/>/{sub(">",">")}1' resfinder-refs.fna > ../resfinder-refs.fna
cd .. && rm -r resfinder

mkdir card && cd $_
wget --no-check-certificate -qO- https://card.mcmaster.ca/download | grep -Eo "download/0/broadstreet[a-zA-Z0-9./?=_-]*" | sort | uniq | tail -n 1 > card-db-version
cardLink=$(sed  's/^/https:\/\/card.mcmaster.ca\//g' card-db-version)
wget --no-check-certificate -O card-db.tar.gz $cardLink
tar -xvf card-db.*
awk '/>/{sub(">",">")}1' nucleotide_fasta_protein_homolog_model.fasta > ../card-refs.fna
cd .. && rm -r card