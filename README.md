# Metagenomic Sequencing Pipeline

Runs all incoming metagenomic sequencing data through a consistent process,
including cleaning, species assignment, and counting human-infecting viruses.

## Table of contents
- [Working with the output](#working-with-the-output)
- [Adding new data](#adding-new-data)
- [Design](#design)
- [Dependencies](#dependencies)

## Status

As of 2023-04-12, handles the data from nine papers.  Currently only handles
short-read paired-end data.

## Working with the output

### Metadata

The metadata is in `dashboard/metadata_*.json`:

* [`metadata_papers.json`](https://github.com/naobservatory/mgs-pipeline/blob/main/dashboard/metadata_papers.json):
  The papers that this data was collected from.  Includes the nucleic acid type
  (though this should move to being a sample-level attribute), a link to the
  paper, and a list of the paper's bioprojects.

* [`metadata_bioprojects.json`](https://github.com/naobservatory/mgs-pipeline/blob/main/dashboard/metadata_bioprojects.json):
  The bioprojects, with which samples they contain.  Some papers have multiple
  bioprojects, though in the cases I've looked at this doesn't seem to
  represent anything?

* [`metadata_samples.json`](https://github.com/naobservatory/mgs-pipeline/blob/main/dashboard/metadata_samples.json):
  Sample-level metadata.  Multi-level location information and sampling date.

These files are generated by
[`prepare_dashboard_data.sh`](https://github.com/naobservatory/mgs-pipeline/blob/main/dashboard/prepare-dashboard-data.sh).

All NCBI taxonomic IDs are relative to the 2022-12-01 release
([zip](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-12-01.zip)).
You can get the scientific name for a taxid with:

    $ cat names.dmp | awk -F'\t' '$1==12234&&$7=="scientific name"{print $3}'
    Tobamovirus

The parent taxid for a taxid with:

    $ cat nodes.dmp | awk -F'\t' '$1==12234{print $3}'
    675071

And the child taxids for a taxid with:

    $ cat nodes.dmp | awk -F'\t' '$3==12234{print $1}'
    12235
    12238
    12239
    12240
    12241
    ...

Or, with names,

    $ for taxid in $(cat nodes.dmp | awk -F'\t' '$3==12234{print $1}'); do
        echo $taxid $(cat names.dmp | awk -F'\t' '$1=='$taxid'&&$7=="scientific name"{print $3}')
      done
    12235 Cucumber green mottle mosaic virus
    12238 Odontoglossum ringspot virus
    12239 Pepper mild mottle virus
    12240 Sunn-hemp mosaic virus
    12241 Tobacco mild green mosaic virus
    ...

### Data

The original data comes from the [Sequencing Read
Archive](https://www.ncbi.nlm.nih.gov/sra) (SRA), and bioproject IDs are SRA
bioproject accessions.  For example, "PRJNA729801" for the Rothman 2021 data.
Sample IDs are SRA run accessions, and are globally unique.  For example,
SRR14530724 belongs to PRJNA729801 (Rothman 2021) and (per
[`metadata_samples.json`](https://github.com/naobservatory/mgs-pipeline/blob/main/dashboard/metadata_samples.json)) represents a 2020-09-10 sample from the
`SB` plant.

Output and intermediate steps for each bioproject are in S3, under
`s3://nao-mgs/[bioprojectid]` (or `s3://nao-restricted/...` for private data).
Available files, and their formats:

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
   * Ex: `SRR21452137.alignments.tsv.gz`,
   * Compressed TSV
   * One record for each read that Bowtie2 was able to map back
     to a known human-infecting virus.  Note that we've set the quality score
     to the minimum and you likely want to filter some of these out based on a
     combination of alignment score and trimmed read length.
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

1. Find the bioproject's accession, so we can load the data.  For example
   https://pubmed.ncbi.nlm.nih.gov/34550753/ is `PRJNA729801`.  Right now we
   can only handle short-read paired-end data.
2. Run `./import-accession.sh [accession]`.  Continue in parallel; this doesn't
   have to finish before you get to step #5, it just needs a small head start.
   * If some files hit errors it's fine to re-run; it skips any files that are
     already complete.
3. Make a directory `bioprojects/[accession]/metadata` and:
   b. Create `bioprojects/[accession]/metadata/name.txt` with the short name of
      the associated paper.  For example, `Rothman 2021`.
   c. Create `bioprojects/[accession]/metadata/metadata.tsv` with a list of the
      sample accessions in the first column and anything else you want to record about the sample in the later
      columns.  If you don't have the rest of the metadata sorted out yet and
      just want to unblock the pipeline you can put the accessions only with:

          curl -sS 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA648796&result=read_run&limit=0' | \
            grep -v ^run_accession | \
            awk '{print $1}' > metadata.tsv

4. Make a directory `papers/AuthorYear/` (matching `name.txt` but without the
   space) and:
   a. Put the paper link in `papers/AuthorYear/link.txt`
   b. Put "RNA" or "DNA" in `papers/AuthorYear/na_type.txt`
5. `./run.py --bioproject=[accession]`
6. Collect full metadata, update metadata.tsv, and modify
   `dashboard/sample_metadata_classifier.py` to handle it.
7. If you leave out any samples, perhaps because they represent something we're
   not interested in, document which samples are included in
   `papers/AuthorYear/subset.txt`.
8. Run `dashboard/prepare_dashboard_data.sh`.

When collecting metadata, put it in this repo under
`bioprojects/[accession]/metadata/`.  Include both the metadata files and the
scripts that prepare it.

Each bioproject has a `bioprojects/[accession]/metadata/metadata.tsv` where the
first column is the sample ID and the remaining columns are
bioproject-specific.  The interpretation of those columns goes in
`dashboard/prepare-dashboard-data.py`.

## Design

Each stage reads from an S3 directory under `s3://nao-mgs/[bioprojectId]/` and
writes to a different one.

### Species Classification

Kraken to assign taxonomic identifiers to reads.  We may want to integrate
Bracken to estimate per-species relative abundance.  We also may need to build
our own database or extend an existing one to ensure we have the level of
coverage we need for human viruses.

When we do want to build our own database kraken2-build has good [documentation
here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases).

## Dependencies

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

The adapter removal step downloads the fastq files one sample at a time from S3 to the local machine.
The step also generates local output files there before copying them to S3.
On Linux these files are stored in `/tmp`.
If you don’t have enough space available in `/tmp`, AdapterRemoval will crash.
To check the available space in `/tmp` run `df -H`.
You should have at least 2x the size of your largest pair of fastq files.

To resize `/tmp`, edit `/etc/fstab`.
If there is an entry for `/tmp`, add the option `size=64G` (or whatever size you need) to the 4th column.
If not, add this line to the end of the file (tab-separated):

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
As of 2023/09/26, there is a [bug](https://github.com/microsoft/onnxruntime/issues/17631) in the latest version of one of the dependencies of RiboDetector,
which causes the program to crash. To avoid this, install an earlier version of the dependency:

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

For now I'm using the Standard pre-built database (see [documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)), with a 16GB cap.  Longer term what
we want depends on how much memory we're ok using; full Standard doesn't quite
fit on a machine with 64GB of RAM.

```
mkdir ~/kraken-db
cd ~/kraken-db
# avoids needing 2x the storage
aws s3 cp s3://genome-idx/kraken/k2_standard_16gb_20221209.tar.gz - | tar -xvz
```

Don't update the version of Kraken's DB without talking to everyone: we would
need to reprocess all older data to handle taxonomy changes and to keep
everything consistent.  Also keep this in sync with the version of the taxonomy
in prepare-dashboard-data.sh.

### BowTie2

#### Install

```
wget -O bowtie2.zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.2/bowtie2-2.5.2-linux-x86_64.zip/download" 
unzip bowtie2.zip 
rm bowtie2.zip
```

#### Set up Bowtie2 database

To create a Bowtie2 database, we need to download genomes from NCBI, using ncbi-genome-download. To run the affiliated script `gimme_taxa.py`, you will also need to install the dependencies `ete3` and `six`. 

```
pip install ncbi-genome-download
python -m pip install ete3 six ncbi-genome-download
```

Now you can run `build_bowtie_db.py`, which will create the bowtie2 database. This will take quite some time so it's best to run this command within a `screen` session to not inadvertantly end the script when closing your terminal.


## Operations

We normally run this pipeline on EC2 instances, generally c6a.8xlarge ones.

### Rerunning for all Bioprojects

The pipeline runs serially within each bioproject. It would be possible to
extend it to be much more parallel, but we haven't felt the need yet.

When we parallize it we usually do it over bioprojects, by running
./reprocess-bioprojects.py and specifying a maximum number of parallel jobs.
For example:

```
mgs-pipeline $ ./reprocess-bioprojects.py 12 prefix --run-arguments-here
```

This will run the pipeline once for each bioproject, including restricted
bioprojects if available in ../mgs-restricted.  If the stages you're running
require a lot of memory or disk, use a number lower than 12, potentially 1.

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
 --bioproject PRJEB49260 --stages hvreads
```

This prints, for every currently executing instance of the pipeline, what
arguments it was started with.

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

4. Follow the instructions above to rerun across all bioprojects.
