import os
import argparse
import tempfile
import subprocess
import shutil
import gzip

# Constants
S3_BUCKET = "s3://nao-mgs"
THISDIR = os.path.abspath(os.path.dirname(__file__))
OUTPUT_BUCKET = "s3://lenni-analysis/ribodetector"


def create_temp_dir(parent_path=THISDIR):
    """Create and return a temporary directory path."""
    return tempfile.mkdtemp(dir=parent_path)


def copy_from_s3_to_local(s3_uri, local_path):
    """Copy file from S3 to local."""
    try:
        subprocess.check_call(["aws", "s3", "cp", s3_uri, local_path], stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print(f"Error copying {s3_uri} to {local_path}")


def copy_inputs_from_s3(bioproject, sample, dest_path):
    """Copy input files from S3 bucket to local temp directory."""
    extensions = [".collapsed.gz", ".collapsed.truncated.gz", ".pair1.truncated.gz", ".pair2.truncated.gz"]
    sample_uris = [os.path.join(S3_BUCKET, bioproject, "cleaned", sample + ext) for ext in extensions]

    for s3_uri in sample_uris:
        local_file = os.path.join(dest_path, os.path.basename(s3_uri))
        copy_from_s3_to_local(s3_uri, local_file)


def calculate_average_read_length(file_path):
    """Calculate and return the average read length of a FASTQ file."""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    total_length = sum(len(line.strip()) for i, line in enumerate(lines) if i % 4 == 1)
    total_reads = len(lines) // 4

    avg_length = int(total_length / total_reads)
    print(f"Average read length of {file_path} is {avg_length:.2f}")
    return avg_length


def count_reads_in_fastq(fq_file):
    """Count and return the number of reads in a FASTQ file."""
    with open(fq_file, "r") as f:
        return sum(1 for _ in f) // 4


def decompress_and_rename(gzip_file):
    """Decompress a gzip file and return fastq filenames."""
    subprocess.check_call(["gunzip", gzip_file])
    unzipped_fname = gzip_file[:-3]
    fastq_fname = gzip_file[:-3] + ".fq"
    os.rename(unzipped_fname, fastq_fname)
    return fastq_fname

def process_single_read(bioproject, sample, read_path, threads, chunk_size):
    """Process a single FASTQ file with ribodetector and report results."""
    fastq_file = decompress_and_rename(read_path + ".gz")
    avg_read_length = calculate_average_read_length(fastq_file)

    # Define Ribodetector output path
    output_path = read_path + ".nonrrna.fq"
    print(f"Running ribodetector_cpu on: {fastq_file}")

    subprocess.check_call([
        "ribodetector_cpu",
        "-t", str(threads),
        "-l", str(avg_read_length),
        "-i", fastq_file,
        "-e", "rrna",
        "--chunk_size", str(chunk_size),
        "-o", output_path
    ])

    report_and_upload_to_s3(fastq_file, output_path, bioproject, sample)

def process_paired_reads(bioproject, sample, read_paths, threads, chunk_size):
    """Process paired FASTQ files with ribodetector."""
    fastq_files = [decompress_and_rename(read + ".gz") for read in read_paths]
    avg_read_length = sum([calculate_average_read_length(f) for f in fastq_files]) // 2

    # Define Ribodetector output path for paired reads
    output_paths = [read_paths[0] + ".nonrrna.fq", read_paths[1] + ".nonrrna.fq"]

    subprocess.check_call([
        "ribodetector_cpu",
        "-t", str(threads),
        "-l", str(avg_read_length),
        "-i", fastq_files[0], fastq_files[1],
        "-e", "rrna",
        "--chunk_size", str(chunk_size),
        "-o", output_paths[0], output_paths[1]
    ])

    report_and_upload_to_s3(fastq_files[0], output_paths[0], bioproject, sample)
    report_and_upload_to_s3(fastq_files[1], output_paths[1], bioproject, sample)


def report_and_upload_to_s3(original_file, processed_file, bioproject, sample):
    """Report results and upload processed file to S3."""
    original_count = count_reads_in_fastq(original_file)
    processed_count = count_reads_in_fastq(processed_file)
    rrna_count = original_count - processed_count
    percentage_rrna = (rrna_count / original_count) * 100

    print(f"Original number of reads in {os.path.basename(original_file)}: {original_count}")
    print(f"Number of reads after ribodetector: {processed_count}")
    print(f"Difference (number of reads classified as rRNA): {rrna_count} "
          f"({percentage_rrna:.2f}% classified as rRNA)")

    # Compress the file using gzip
    compressed_file = processed_file[:-3]+ ".gz"
    with open(processed_file, 'rb') as orig_file, gzip.open(compressed_file, 'wb') as zipped_file:
        shutil.copyfileobj(orig_file, zipped_file)

    # Upload the compressed results to S3
    s3_path = os.path.join(OUTPUT_BUCKET, bioproject, sample, os.path.basename(compressed_file))
    subprocess.check_call(["aws", "s3", "cp", compressed_file, s3_path], stdout=subprocess.DEVNULL)
    

def run_ribodetector_on_sample(bioproject, sample, temp_path, threads=20, chunk_size=256):
    """Run ribodetector on given sample's reads."""
    single_read_extensions = [".collapsed", ".collapsed.truncated"]
    paired_read_extensions = [".pair1.truncated", ".pair2.truncated"]

    single_reads = [os.path.join(temp_path, sample + ext) for ext in single_read_extensions]
    paired_reads = [os.path.join(temp_path, sample + ext) for ext in paired_read_extensions]

    for read in single_reads:
        process_single_read(bioproject, sample, read, threads, chunk_size)
    process_paired_reads(bioproject, sample, paired_reads, threads, chunk_size)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run RiboDetector')
    parser.add_argument('--bioproject', required=True, help='Bioproject ID to process.')
    parser.add_argument('--sample', required=True, help='SRA run accession of sample to process.')
    return parser.parse_args()


def main():
    args = parse_arguments()
    temp_directory = create_temp_dir()
    try:
        copy_inputs_from_s3(args.bioproject, args.sample, temp_directory)
        run_ribodetector_on_sample(args.bioproject, args.sample, temp_directory)
    finally:
        shutil.rmtree(temp_directory)


if __name__ == "__main__":
    main()
