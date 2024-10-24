# Metagenomic Sequencing Pipeline

Runs all incoming metagenomic sequencing data through a consistent process,
including cleaning, species assignment, and counting human-infecting viruses.

## Table of contents
- [Working with the output](#working-with-the-output)
- [Adding new data](#adding-new-data)
- [Design](#design)
- [Dependencies](#dependencies)

## Status

Deprecated

### Data

[deprecated]

Output and intermediate steps for each delivery are in S3, under
`s3://nao-mgs/[deliveryid]` (or `s3://nao-restricted/...` for private data).
Available files, and their formats:

(Note that directories downstream from kraken have a -YYYY-MM suffix,
representing the version of the kraken DB they were generated with.)

* `raw/`: Input, as downloaded from SRA
   * ex: `SRR14530724_1.fastq.gz` and `SRR14530724_2.fastq.gz`
   * Gzipped FastQ format
   * Paired-end reads split across two parallel files, one for the forward
     reads and another of exactly the same length for the reverse reads.

* `cleaned/`: Intermediate, output of `AdapterRemoval2`
   * ex: SRR14530724.collapsed.gz
   * Gzipped FastQ format
   * Multiple files, depending on what `AdapterRemoval2` did
     * `collapsed` if it was able to combine the pair into one, `pair1` and
       `pair2` otherwise.
     * `truncated` if it did any quality trimming.
     * `singleton` if it dropped either the forward or reverse read.
     * `discarded` if it dropped the whole pair
     * Plus a `settings` file with info on how cleaning went

* `ribofrac/`: Intermediate, fraction of rRNA from fast `RiboDetector`
    * ex: COMO1.ribofrac.txt
    * txt
    * Fraction of rRNA in sample using subset of input reads

* `riboreads/`: Intermediate, read IDs of `RiboDetector` output
    * ex: SRR13167436.riboreads.txt.gz
    * txt
    * Stores rRNA read IDs

* `processed/`: Intermediate, output of `Kraken2`
   * ex: `SRR14530724.collapsed.kraken2.tsv.gz`
   * Gzipped TSV
   * Kraken output format
     * One line per input read with taxonomic classification and detailed hit
       information.

* `samplereads/`: Itermediate, example reads by taxonomic category
  * ex: `SRR14530724.sr.tsv.gz`
  * TSV
  * 100,000 examples (or as many as available) from each sample of:
    * a: all reads
    * b: bacterial reads
    * v: viral reads
    * h: human viral reads
  * These are sampled by assuming that order within a sequencing run doesn't
    matter: it takes the first 100,000 of each category.
  * Columns are:
    * Category (a/b/v/h)
    * Read ID

* `readlengths`: Output, distribution of read lengths by taxonomic category
  * ex: `SRR14530724.rl.json.gz`
  * JSON
  * category -> length -> count
  * Categories:
    * a: all reads
    * b: bacterial reads
    * v: viral reads
    * h: human viral reads
  * Length is a number, and then "NC" is cases where cleaning was unable to
    collapse the read.  In these cases all we know is that the fragment was
    long enough to not get overlap.
  * Because we currently don't have a stage that does adapter removal without
    quality trimming, these will be an underestimate of fragment lengths,
    especially for lower quality reads.

* `allmatches/`: Intermediate, subset of `Kraken2` output matching human viruses
   * ex: `SRR14530724.allmatches.tsv`
   * TSV
   * One row for each read with any hits against a human-infecting virus.
   * Kraken output format
   * Match can be in the assigned taxid, or any of the hits.

* `hvreads/`: Output, data that backs the [read
  viewer](https://www.jefftk.com/mgs-counts/reads)
   * Ex: `SRR14530724.hvreads.json`
   * One record for each record in allmatches, joined back to the cleaned
     sequence and quality data.
   * JSON
   * Read ID to Kraken output and cleaned read

* `alignments/`: Output, alignment data that will later back the dashboard.
   * Ex: `SRR21452137.hv.alignments.tsv.gz`,
         `SRR21452137.human.alignments.tsv.gz`
   * Compressed TSV
   * One record for each read that Bowtie2 was able to map back to a genome in
     its DB.
   * We run this twice: once with a standard human genome DB ("human"), and
     again with a custom DB of human-infecting viruses ("hv").
     * For Human reads it runs with default setting, and the score is zero for
       perfect matches and increasingly negative for worse matches.
     * For HV reads we're running with custom settings, and the score is
       positive, higher for better matches.  Some reads with positive scores
       are quite low quality matches: you likely want to filter with a
       combination of alignment score and trimmed read length (ex:
       Simon's been doing `score/ln(length) > 22`).
   * Non-collapsed reads will appear twice, one for the forward read and then
     one for the reverse.
   * Columns:
     * Read ID
     * Best-match genome
     * Best-match taxid
     * CIGAR string
     * Genome start position
     * Alignment score
     * Trimmed read length

* `humanviruses/`: Output, data that backs the
  [dashboard](https://www.jefftk.com/mgs-counts/).
   * Ex: `SRR14530724.humanviruses.tsv`
   * TSV
   * One row for each human-infecting virus observed in the sample
   * Columns are:
     * Taxid: NCBI Taxonomic ID, from the 2022-12-01 taxdmp release
     * Count: how many reads Kraken assigned to this taxid
     * Scientific name: the NCBI scientific name for this taxid

* `cladecounts/`: Output, taxonomic counts
   * ex: `SRR21452136.tsv.gz`
   * Gzipped TSV
   * Collated Kraken output
   * Columns:
     * Taxid: NCBI Taxonomic ID, from the 2022-12-01 taxdmp release
     * Direct assignments: how many reads in this sample Kraken assigned to
       this taxid.
     * Direct hits: how many reads in this sample had any 35-mer matches to the
       Kraken DB for this taxid
     * Clade assigments: how many reads in this sample Kraken assigned to this
       taxid or anywhere in its clade.
     * Clade hits: how many reads in this sample had any 35-mer matches to the
       Kraken DB for this taxid or anywhere in its clade.
   * For example, Tomaboviruses (12234) in `SRR21452136.tsv.gz` is:
         12234    920   45914   954843   1027715
     This is saying that:
     * 920 reads were assigned to "Tobamovirus".
     * 45,914 reads had some 35-mer that was common to Tobamoviruses and
       nothing more specific.
     * 954,843 reads were assigned to Tobamovirus or something within it's
       clade (ex: PMMoV).
     * 1,027,715 reads had at least one 35-mer hit within this clade.

## Adding new data

We are no longer importing new public data.

## Design

Each stage reads from an S3 directory under `s3://nao-mgs/[deliveryId]/` and
writes to a different one.

### Species Classification

Kraken to assign taxonomic identifiers to reads.

## Dependencies

### E-utilities (for working with NCBI)
See the NCBI [documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/). 
```
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}
sudo yum install perl-Time-HiRes

```

### AdapterRemoval2

Visit https://github.com/MikkelSchubert/adapterremoval/ to get the latest
version number, and then:

```
wget -O adapterremoval-2.3.3.tar.gz \
     https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.3.tar.gz
tar xvzf adapterremoval-2.3.3.tar.gz
cd adapterremoval-2.3.3
make
sudo make install

```

The adapter removal step downloads the fastq files one sample at a time from S3
to the local machine.  The step also generates local output files there before
copying them to S3.  On Linux these files are stored in `/tmp`.  If you don’t
have enough space available in `/tmp`, AdapterRemoval will crash.  To check the
available space in `/tmp` run `df -H`.  You should have at least 2x the size of
your largest pair of fastq files, multiplied by the maximum number of files you
will need to process simultaneously.

To resize `/tmp`, edit `/etc/fstab`. If there is an entry for `/tmp`, add the
option `size=64G` (or whatever size you need) to the 4th column.  If not, add
this line to the end of the file (tab-separated):

```
tmpfs  /tmp  tmpfs  size=64G  0  0
```

then run:

```
sudo mount -o remount /tmp/
```

### RiboDetector
See [documentation](https://github.com/hzi-bifo/RiboDetector). To install:

```
pip install ribodetector
```

As of 2023-08-26, there is a
[bug](https://github.com/microsoft/onnxruntime/issues/17631) in the latest
version of one of the dependencies of RiboDetector, which causes the program to
crash. To avoid this, install an earlier version of the dependency:

```
pip install onnxruntime==1.15.1
```

### Kraken

#### Install

```
git clone git@github.com:DerrickWood/kraken2.git
cd kraken2/
./install_kraken2.sh ~/kraken2-install
```

#### Set up Kraken database

We're using the Standard pre-built database (see
[documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)).
The prepare-shm-kraken.sh script automatically downloads and sets it up, but
first you need to configure RAM.  Edit `/etc/fstab` and add:

```
none     /dev/shm      tmpfs  defaults,size=85G      0 0
```

The 85G comes from starting with the 78GB listed for the Standard kraken DB at
https://benlangmead.github.io/aws-indexes/k2 and then adding 6.5GB for the
custom human virus and human bowtie DBs.

See "Updating the Taxnonomy and Kraken DB" below if it's out of date and you'd like
it not to be.

### BowTie2

#### Install

```
wget -O bowtie2.zip "https://sourceforge.net/deliverys/bowtie-bio/files/bowtie2/2.5.2/bowtie2-2.5.2-linux-x86_64.zip/download" 
unzip bowtie2.zip 
rm bowtie2.zip
```

#### Download pre-built human genome database

For detecting human reads we use the standard pre-built telomere-to-telomere
"Human / CHM13plusY" database from
https://benlangmead.github.io/aws-indexes/bowtie.  See
https://www.science.org/doi/10.1126/science.abj6987 for the construction of
this genome.

```
cd mgs-pipeline/bowtie
aws s3 cp s3://genome-idx/bt/chm13.draft_v1.0_plusY.zip .
unzip chm13.draft_v1.0_plusY.zip
mv chm13.draft_v1.0_plusY/* .
rmdir chm13.draft_v1.0_plusY
rm chm13.draft_v1.0_plusY.zip
```

#### Build custom human viral database

To create a Bowtie2 database, we need to download genomes from NCBI, using
ncbi-genome-download. To run the affiliated script `gimme_taxa.py`, you will
also need to install the dependencies `ete3` and `six`.

```
pip install ncbi-genome-download
python -m pip install ete3 six ncbi-genome-download
```

Now you can run `build_bowtie_db.py`, which will create the bowtie2
database. This will take quite some time so it's best to run this command
within a `screen` session to not inadvertantly end the script when closing your
terminal.


## Operations

We normally run this pipeline on c6a.16xlarge EC2 instances.  We need 128GB+ of
memory; the binding constraint is that the Standard Kraken DB is 78GB.

### Rerunning for all Deliveries

While ./run.py runs serially for a single delivery, we can run in parallel at
the sample level with ./reprocess.py. You specify a maximum number
of parallel jobs.  For example:

```
mgs-pipeline $ ./reprocess.py \
    --max-jobs 12 \
    --log-prefix prefix \
    --sample-level \
    --shuffle
    -- \
    --run-arguments-here \
```

This will run the pipeline once for each delivery, including restricted
deliveries if available in ../mgs-restricted.  If the stages you're running
require a lot of disk, use a number lower than 12, potentially 1: the script
doesn't know how "heavy" a job is, and will run out of disk space if you tell
it to do do too much in parallel.

Job output is under log/ in files named by the date and the prefix you supply.
So if I ran the above on 2023-01-01 I'd expect to see files like:

```
...
log/2023-01-01.prefix.PRJNA924011
log/2023-01-01.prefix.PRJNA943189
log/2023-01-01.prefix.PRJNA966185
...
```

If the job fails, the last line in the log file will start with "ERROR:" and
then have the exit code.

### Screen Oversight

You can check in on parallelized jobs under screen with:

```
mgs-pipeline $ pipeline-operation/screen-summary.py
8769..assembly:
4.collapsed.gz
hvreads: handling ERR7850094.collapsed.truncated.gz in /tmp/tmphu0ym26c
download: s3://nao-mgs/PRJEB49260/cleaned/ERR7850094.collapsed.truncated.gz to
.
/ERR7850094.collapsed.truncated.gz
hvreads: handling ERR7850094.discarded.gz in /tmp/tmpegmh_drc
download: s3://nao-mgs/PRJEB49260/cleaned/ERR7850094.discarded.gz to
./ERR785009
4.discarded.gz
hvreads: handling ERR7850094.pair1.truncated.gz in /tmp/tmpmel9okco
download: s3://nao-mgs/PRJEB49260/cleaned/ERR7850094.pair1.truncated.gz to
./ERR
7850094.pair1.truncated.gz
````

This prints the bottom ten lines of each active screen session on the computer.

You can also run:

```
mgs-pipeline $ pipeline-operation/print-running-jobs.sh
 --delivery PRJEB49260 --stages hvreads
```

For every currently executing instance of ./run.py this prints the arguments it
was started with.

### Regenerating Data

Normally the pipeline doesn't repeat work, but when the code changes some data
typically also needs to change.

Normally the flow is:

1. Run the pipeline on a single sample (`--sample RUN_ACCESSION`) until you're
   happy with what it does.

2. Find where it checks whether your output already exists. For example, for
   `hvreads` there's a line like:
   `existing_outputs = get_files(args, "hvreads")`.

3. Add a `min_date` argument to the `get_files` call, like
   `min_date='2023-09-18'`.

4. Follow the instructions above to rerun across all deliveries.

### Updating the Taxnonomy and Kraken DB

Before starting, check in with all other users of the pipeline, and ensure no
one wants you to hold off.  Also check how recent the databases are on
https://benlangmead.github.io/aws-indexes/k2.  If they're not recent, it may be
worth waiting for the next update.

Once you've decided to go ahead:

1. Move human-viruses-raw.tsv, human-viruses.tsv, plus, within dashboard,
   *.dmp, top_species_counts, and top_species_scratch into a scratch
   location. Delete hvreads, readlengths, ribofrac, cladecounts and
   allmatches.

2. Update download-taxonomy.sh to pull the latest taxonomy.  Run the script.

3. Run download-human-viruses.sh

4. Put the year and month in reference-suffix.txt

5. Reprocess the deliveries we need to have up to date.  The command will look
   something like:

```
./reprocess.py \
   --log-prefix kraken-update \
   --max-jobs 8 \
   --shuffle \
   --sample-level \
   --deliveries DELIVERY_1,DELIVERY1_,...DELIVERY_N \
   --
```

6. Manually compare the output to the previous generation.

7. Consider whether we need to keep the previous generation data (or previous
   previous etc) around.
