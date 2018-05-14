#install local copy of samtools
wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
tar -xvjf samtools-1.8.tar.bz2
cd samtools-1.8
make
cd ..