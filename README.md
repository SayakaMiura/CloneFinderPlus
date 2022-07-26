# CloneFinder+_v0.1

## Description
CloneFinder+ is a method aimed at estimating individual clone genotypes and frequencies within a tumor sample using a phylogenetic approach. See Huzar et al. (ref. 1) for the detail. The CloneFinder+ program has been developed by Sudhir Kumar. It is written in Python. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below).

## Dependencies
CloneFinder+ is a python script developed in a Windows 64-bit architecture.
1. python 3 (v3.8 was tested)
 > python packages: 
 >>  SciPy
 >>
 >>  NumPy 
 >>
 >>  Biopython 
 >>
 >>  pandas 
 >>
 >>  networkx  
 >>  
 >>  matplotlib 
 >>
 >>  logbook
 >>
 >>  fire 
 >>
 >> pydot
 >>
 >Note: If the installation of these python packages is not easy, you may want to use Anaconda for Python 3 (https://www.anaconda.com/distribution/). Or you can try python -m pip install [package name].

## How to use CloneFinder+

`python Bootstrap_CloneFinderPlus.py [read count file] [tumor site file (optional)] [Number of bootstrap replicates] `

#### input 
- read count file
 
 The input file is a tab-delimited text file, which contains the read counts of wild-type and mutant alleles for each sample. Normal sample should not be included. Each line in the input files gives information for each variant. 
 
 `"XX:ref": Reference read count for the sample, XX
"XX:alt": Variant read count for the sample, XX`
 
- tumor site file (optional)
 
 When migration paths are selected to be inferred, the information of tumor sites need to be provided for each tumor sample. All tumor samples should be listed in the first column, and their site name should be listed in the second column. For example,
 
 `Sample	Site
PrimT	P
Met2T	Met2
…`

#### output files
 - Consensus sequences (fasta file)
 
 Mutant and wild type bases are indicated by “T” and “A,” respectively. 
 - Bootstrap support for each SNV assignment (txt file)
 
For each SNV position, bootstrap support for mutant base assignment is given.

- (optional) Bootstrap support for migration paths (txt file)
 
 For each bootstrap migration path, bootstrap support is given.


#### example
 To run example datasets with 30 bootstrap replicates in Example directory, run the following command.
 
`python Bootstrap_CloneFinderPlus.py Example\ReadCount.txt Example\TumorSite.txt 30`

 The output files will be stored in "Example" directory.

To run without the bootstrap option.
 
`python Bootstrap_CloneFinderPlus.py Example\ReadCount.txt Example\TumorSite.txt 0`

To run example datasets without migration path inference.
 
`python Bootstrap_CloneFinderPlus.py Example\ReadCount.txt NA 30`

`python Bootstrap_CloneFinderPlus.py Example\ReadCount.txt NA 0`


## Reference:
[1] Jared Huzar, Madelyn Shenoy, Maxwell Sanderford, and Sudhir Kumar, and Sayaka Miura. Reliable tumor evolution  estimates using bulk sequencing data (2022) In preparation


## Copyright 2022, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
