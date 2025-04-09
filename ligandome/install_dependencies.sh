#!/bin/bash

WORKING_DIR=$(pwd)
mkdir $WORKING_DIR/ligandome/third_party_tools

# install BLAST
if [ "$(uname)" == "Darwin" ]; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-x64-macosx.tar.gz -P ./ligandome/third_party_tools/
    tar zxvpf ./ligandome/third_party_tools/ncbi-blast-2.14.1+-x64-macosx.tar.gz -C ./ligandome/third_party_tools/

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-x64-linux.tar.gz -P ./ligandome/third_party_tools/
    tar zxvpf ./ligandome/third_party_tools/ncbi-blast-2.14.1+-x64-linux.tar.gz -C ./ligandome/third_party_tools/
fi
export PATH=$PATH:$WORKING_DIR/ligandome/third_party_tools/ncbi-blast-2.14.1+/bin

# create BLAST databases
mkdir ./database_exports/blast_dbs/
makeblastdb -in ./ligandome/database_exports/UniProtCanonicalProteome_2023-08-29.fasta -dbtype prot -parse_seqids -out ./ligandome/database_exports/blast_dbs/UniProtCanonicalProteome

# install SACHICA
wget http://research.nii.ac.jp/~uno/code/sachica34.zip  -P ./ligandome/third_party_tools/
mkdir ./ligandome/third_party_tools/SACHICA
unzip ./ligandome/third_party_tools/sachica34.zip -d ./ligandome/third_party_tools/SACHICA
cd ./ligandome/third_party_tools/SACHICA; make; cd $WORKING_DIR
export PATH=$PATH:$WORKING_DIR/ligandome/third_party_tools/SACHICA

# manual install of NetMHCPan 4.1
unzip ./ligandome/third_party_tools/netMHCpan-4.1.zip -d ./ligandome/third_party_tools/
if [ "$(uname)" == "Darwin" ]; then
    sed -i '' -e "s|PATH_TO_NETMHCPAN|$WORKING_DIR/ligandome/third_party_tools|" ./ligandome/third_party_tools/netMHCpan-4.1/netMHCpan

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    sed -Ei "s|PATH_TO_NETMHCPAN|$WORKING_DIR/ligandome/third_party_tools|g" ./ligandome/third_party_tools/netMHCpan-4.1/netMHCpan
fi
export PATH=$PATH:$WORKING_DIR/ligandome/third_party_tools/netMHCpan-4.1