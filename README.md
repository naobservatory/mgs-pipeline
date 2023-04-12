# Metagenomic Sequencing Pipeline

Runs all incoming metagenomic sequencing data through a consistent process,
including cleaning, species assignment, and counting human-infecting viruses.

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
`s3://nao-mgs/[bioprojectid]`.  Available files, and their formats:

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

* `processed/`: Intermediate, output of `Kraken2`
   * ex: `SRR14530724.collapsed.kraken2.tsv.gz`
   * Gzipped TSV
   * Kraken output format
     * One line per input read with taxonomic classification and detailed hit
       information.

* `allmatches/`: Intermediate, subset of `Kraken2` output matching human viruses
   * ex: `SRR14530724.allmatches.tsv`
   * TSV
   * One row for each read with any hits against a human-infecting virus.
   * Kraken output format
   * Match can be in the assigned taxid, or any of the hits.

* `hvreads/`: Output, data that backs the [read
  viewer](https://www.jefftk.com/mgs-counts/reads)
   * Ex: `SRR14530724.hvreads.json`
   * JSON
   * Read ID to Kraken output and cleaned read

* `humanviruses/`: Output, data that backs the
  [dashboard](https://www.jefftk.com/mgs-counts/).
   * Ex: `SRR14530724.humanviruses.tsv`
   * TSV
   * One row for each human-infecting virus observed in the sample
   * Columns are:
     * Taxid: NCBI Taxonomic ID, from the 2022-12-01 taxdmp release
     * Count: how many reads Kraken assigned to this taxid
     * Scientific name: the NCBI scientific name for this taxid

* `allcounts/`: Intermediate, taxonomic counts
   * ex: `SRR14530724.tsv.gz`
   * Gzipped TSV
   * Collated Kraken output
   * Columns:
     * Taxid: NCBI Taxonomic ID, from the 2022-12-01 taxdmp release
     * Assignments: how many reads in this sample Kraken assigned to this taxid
     * Hits: how many reads in this sample had any 35-mer matches to the Kraken
       DB for this taxid
   * For example, Tomaboviruses (12234) is:
         12234      270     105908
     This is saying that only 270 reads were assigned to "Tobamovirus" but
     105,908 reads had some 35-mer that was common to Tobamoviruses.

* `childcounts/`: Output, taxonomic counts including children
   * ex: `SRR14530724.tsv.gz`
   * Gzipped TSV
   * Collated Kraken output, plus columns for children
   * Columns:
     * First three columns are the same as `allcounts/` above.
     * Clade assigments: how many reads in this sample Kraken assigned to this
       taxid or anywhere in its clade.
     * Clade hits: how many reads in this sample had any 35-mer matches to the
       Kraken DB for this taxid or anywhere in its clade.
   * For example, Tomaboviruses (12234) is:
         12234      270     105908  333363  465886
     The second two columns are saying that 333,363 reads were assigned to
     the Tobamovirus clade (ex: assignment to PMMoV) and there were 465,886
     35-mer matches where a read matched anything in the Tobamovirus clade.
     * Hits within the clade isn't a great statistic because there can be
       double counting: one read could have a hit for both Tobamovirus and
       PMMoV.  It would probably be better to count each read toward each clade
       only once.  If this is a statistic we want to work with, let me know and
       I can calculate this from the Kraken output.

## Adding new data

1. Find the bioproject's accession, so we can load the data.  For example
   https://pubmed.ncbi.nlm.nih.gov/34550753/ is `PRJNA729801`.
2. Run `./import-accession.sh [accession]`.  Continue in parallel; this doesn't
   have to finish before you get to step #5, it just needs a small head start.
3. Make a directory `bioprojects/[accession]/metadata` and:
   a. Create `bioprojects/[accession]/metadata/study.json` with contents
      `{"is_paired_end": true}`.  If the bioproject isn't actually paired end then
      put `false`, but you're going to need to do a lot of updating `run.py` to
      handle this case.
   b. Create `bioprojects/[accession]/metadata/name.txt` with the short name of
      the associated paper.  For example, `Rothman 2021`.
   c. Create `bioprojects/[accession]/metadata/metadata.tsv` with a list of the
      sample accessions in the first column and anything else in the later
      columns.  If you don't have the rest of the metadata sorted out yet and
      just want to unblock the pipeline you can put the accessions only:

          aws s3 ls s3://nao-mgs/[accession]/raw/ | \
            awk '{print $NF}' | \
            grep _1 | \
            sed s/_1.fastq.gz// > metadata.tsv
4. Make a directory `papers/AuthorYear/` (matching `name.txt` but without the
   space) and:
   a. Put the paper link in `papers/AuthorYear/link.txt`
   b. Put "RNA" or "DNA" in `papers/AuthorYear/na_type.txt`
5. `./run.py --bioproject=[accession]`
6. Collect full metadata, update metadata.tsv, and modify
   `dashboard/prepare-dashboard-data.py` to handle it.

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
wget -O adapterremoval-2.3.1.tar.gz \
     https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.1.tar.gz
tar xvzf adapterremoval-2.3.1.tar.gz
cd adapterremoval-2.3.1
make
sudo make install
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
