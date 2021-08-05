# Rmapper

Rmapper is a tool for optical map assembly using de Bruijn graph.

## Compilation

You can compile `rmapper` executing `make` in the main directory.

## Usage
After compiling, you can build the bi-labelled de Bruijn graph using

``` bash
./rmapper-build infile [-k ksize] [-d dsize] [-f tf] [-l tl] [-r rev] [-o output]

Computes the bi-labelled de Bruijn graph of the error-corrected rmaps in [infile].
  ksize: [integer] - k-mer size (def. 6)
  dsize: [integer] - minimum distance between two k-mers in a bi-label. (def. 15000)
     tf: [integer] - fragment proximal error tollerance. (def. 500)
     tl: [integer] - length proximal error tollerance. (def. 2000)
    rev: [boolean] - include reverse rmaps. (def. false)
 output: [string]  - output directory that must exists. (def. .)
```
Then the graph can be traversed using

``` bash
./rmapper-traverse indir

Traverse the bi-labelled de Bruijn graph in the directory [indir].
```

## Input format

The input format `[NAME_FILE].val` stores each rmap in three lines. The first line contains the Rmap ID, the second line contains the restriction enzyme used, followed by the Rmap fragments in kbp.

Example:

``` 
Rmap_0
        27 enzyme 3.262 7.618 11.667 15.579 10.687 7.727 26.177 3.523 7.565 18.481 4.317 28.016 11.942 21.72 27.698 8.206 1.521 4.169 4.631 4.052 10.489 10.098 3.471 6.938 3.137 16.248 9.69 

Rmap_1
        24 enzyme 4.476 2.264 16.449 16.598 1.881 7.788 12.494 16.856 4.746 3.001 12.706 2.923 6.043 1.639 16.368 7.555 14.451 6.346 9.003 9.625 10.428 7.868 5.059 8.693 

Rmap_2
        24 enzyme 11.047 20.374 1.801 9.055 4.151 6.723 6.537 4.25 4.446 3.899 4.517 9.679 16.132 3.104 6.893 5.008 9.827 10.458 4.102 4.647 4.084 1.537 8.192 9.935 

```

## Examples

You can download the example data from [here](https://drive.google.com/drive/folders/1CxYtVUTEYq4CybTHCxfyj9E4jNCWd_Xt?usp=sharing) using

### *E. coli* simulated rmaps
Download the rmaps
``` bash
wget -O ecoli_corrected.val 'https://docs.google.com/uc?export=download&id=1SUt7plSAZ6d7uHofxnIk1qrYUWbAzQJt'
```
Build and traverse the graph
``` bash
mkdir ecoli_output
./rmapper-build ecoli_corrected.val -k 6 -d 15000 -f 500 -o ecoli_output/
./rmapper-traverse ecoli_output/
```

<!--- 
wget -O ecoli_sim.out 'https://docs.google.com/uc?export=download&id=1QcU_3R4m3YQaTQ5Vj7ByTQjjcpj1dV-K'
wget -O ecoli.bnx 'https://docs.google.com/uc?export=download&id=1Erd0WlRnHhvtvAkk4Ly-7oHmI_iEzH9-' 
--->

### *Human* simulated rmaps 
Download the rmaps
``` bash
wget -O human_corrected.val 'https://docs.google.com/uc?export=download&id=1mCXT67lwB2Zh0FCnoIX_DTP2k7ZK51vE'
```
Build and traverse the graph
``` bash
mkdir human_output
./rmapper-build ecoli_corrected.val -k 6 -d 25000 -f 1500 -o human_output/
./rmapper-traverse human_output/
```
<!---
wget -O human_sim.out 'https://docs.google.com/uc?export=download&id=1zUTJr6hv3Lj5deUG2tmWJoRL-zkCfpuH'
wget -O human.bnx 'https://docs.google.com/uc?export=download&id=1Hvr3vnhTdTd_GENg81XEw5ID5G9MYpWL'
--->


## Citation 

Please, if you use this tool in an academic setting cite the following paper:

    @inproceedings{Mukherjee20,
    author    = {Kingshuk Mukherjee and 
                 Massimiliano Rossi and 
                 Leena Salmela and 
                 Christina Boucher},
    title     = {Fast and efficient Rmap assembly using bi-labelled de Bruijn graph},
    booktitle = {20th International Workshop on Algorithms in Bioinformatics (WABI 2020)},
    year      = {2020},
    series    = {LIPIcs},
    volume    = {172},
    pages     = {9:1--9:16},
    publisher ={Schloss Dagstuhl-Leibniz-Zentrum f{\"u}r Informatik},
    url       = {https://doi.org/10.4230/LIPIcs.WABI.2020.9}
    }

## Authors

### Theoretical results:

* Kingshuk Mukherjee
* Massimiliano Rossi
* Leena Salmela
* Christina Boucher

### Implementation:

* [Kingshuk Mukherjee](https://github.com/kingufl)

### Documentation:

* [Massimiliano Rossi](https://github.com/maxrossi91)