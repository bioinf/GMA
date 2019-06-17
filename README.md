# Genome Mappability Score Analyzer

Project Source Repository: [https://sourceforge.net/projects/gma-bio/](https://sourceforge.net/projects/gma-bio/)



## Quick Start

### Installation

```bash
git clone https://github.com/bioinf/GMA.git && cd GMA
chmod +x bin/* wrapper.sh recordmem.sh
```

### Usage

```bash
./wrapper.sh [fasta file] [Optional arguments]
```

### Optional arguments

```
-l   Read length. Default: 100
-s   Error rate. Default: 0.00
-o   Insert size (for paired end reads). Default: 0
```

### Results

`results_m.d.Y_H.M.S/mapred.txt` — Mappability track.  
`results_m.d.Y_H.M.S/memory_used.log` — Every second dump of the amount of all memory allocated for a given user (during script execution)
