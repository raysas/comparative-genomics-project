#!/bin/sh

# ------------- space to downlaod tools in ------------------
mkdir -p test/tools/
cd test/tools/

# ------------ pipeline 1: BLAST+ and MCL installation ----------------
sudo apt install ncbi-blast+
sudo apt-get install mcl           

# ------------ MCScanX installation ----------------
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX
make
cd ..

# ------------ circos installation ----------------
curl -O https://circos.ca/distribution/circos-0.69-10.tgz
tar xvfz circos-0.69-10.tgz 
cd circos-0.69-10
# -- need to check if perl is installed
which perl
./circos -modules # need to install all missing modules
# -- making it global
# echo 'export PATH="/mnt/c/Users/rayan/Documents/saclay/courses/comparative-genomics/project/test/tools/circos-0.69-10/bin:$PATH"' >> ~/.bashrc
# source ~/.bashrc
cd ..

# ------------ MUMmer installation ----------------
# wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
# tar -xzvf mummer-4.0.0rc1.tar.gz
# cd mummer-4.0.0rc1
# ./configure
# make
# sudo make install
# cd ..
conda remove mummer -y
conda install -c bioconda mummer
sudo apt install gnuplot

# ------------- R packages ------------------
conda install --solver=classic -c conda-forge r-ggplot2 r-dplyr r-rcolorbrewerg