[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://img.shields.io/packagist/l/doctrine/orm.svg)
# Bioinformatics project
[Bioinformatics course](http://www.fer.unizg.hr/predmet/bio) project: Improving Bloom Filter Performance on Sequence Data Using k-mer Bloom Filters
## Build dependencies
Cmake is needed to build this project. Install as follows: 
```
sudo apt-get install cmake
```
## Build instructions 
To build this project run following commands: 
```
git clone git@github.com:Daria2002/bioinformatics_project.git
cd bioinformatics_project
mkdir build
cd build
cmake ..
make
```
This will automatically download, build and instal the required Bloom filter [library](https://github.com/mavam/libbf/).
## Run instructions 
Download any fasta file from the following [link](http://bacteria.ensembl.org/index.html).
Run the executable from the build folder as follows:
```
./main file.fasta kmerLength outputFile.txt
```
or with the default output file name results.txt in results folder:
```
./main file.fasta kmerLength
```
