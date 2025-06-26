import sys
import argparse
import logging
import pandas as pd
from Bio import SeqIO
from glob import glob


def setup_logging(logfile=None, level=logging.INFO):
    logger = logging.getLogger()
    logger.setLevel(level)

    # Clear any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    handler = logging.FileHandler(
        logfile, mode='w') if logfile else logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter('[%(levelname)s] %(message)s')
    handler.setFormatter(formatter)

    logger.addHandler(handler)


def reformat_fasta_headers(
        accession_table, taxid_mapping, database_dir,
        output_fasta, output_report, log=None):
    setup_logging(logfile=log)
    logging.info("Starting to reformat fasta headers.")

    # Load accession-taxid table
    df = pd.read_csv(accession_table, sep="\t", header=None,
                     names=['accession'], index_col=0)
    mapping = pd.read_csv(taxid_mapping, sep="\t", index_col='accession')
    df = df.join(mapping, how='left')
    logging.info(f"Loaded {len(df)} accessions with taxids.")
    # Prepare output
    report_rows = []
    seq_counter = 1
    output_records = []
    logging.info(f"Output fasta file: {output_fasta}")
    logging.info(f"Output report file: {output_report}")
    logging.info(f"Database directory: {database_dir}")

    with open(output_fasta, "w") as f_out:
        for acc, row in df.iterrows():
            taxid = row['taxid']
            fna_path = glob(f"{database_dir}/{acc}/{acc}*.fna")
            if not len(fna_path):
                logging.warning(f"Warning: {acc} fasta not found, skipping.")
                continue
            for fn in fna_path:
                # Parse and rewrite headers
                for record in SeqIO.parse(fn, "fasta"):
                    new_id = f"{seq_counter}-{taxid}"
                    report_rows.append(
                        {
                            'integer': seq_counter,
                            'taxid': taxid,
                            'original_header': record.description})
                    record.id = new_id
                    record.description = ""
                    SeqIO.write(record, f_out, "fasta")
                    seq_counter += 1

    # Write metadata report
    pd.DataFrame(report_rows).to_csv(output_report, sep="\t", index=False)
    logging.info(f"Reformatted {seq_counter - 1} sequences.")
    logging.info(f"Fasta file written to {output_fasta}")
    logging.info(f"Metadata report written to {output_report}")
    logging.info("Reformatting complete.")


if __name__ == "__main__":

    try:
        snakemake
        reformat_fasta_headers(
            snakemake.input.accessions,
            snakemake.input.taxid_mapping,
            snakemake.params.directory,
            snakemake.output.fasta,
            snakemake.output.report,
            log=snakemake.log if snakemake.log else None
        )

    except NameError:

        parser = argparse.ArgumentParser(
            description="Reformat fasta headers to include taxid and generate report.")
        parser.add_argument("accessions",
                            help="TSV with accessions")
        parser.add_argument("taxid_mapping",
                            help="TSV mapping accessions to taxid")
        parser.add_argument("--database-dir", required=True,
                            help="Path to base directory containing accession subfolders")
        parser.add_argument("--output-fasta", required=True,
                            help="Path to output fasta file")
        parser.add_argument("--output-report", required=True,
                            help="Path to output report file with metadata")
        parser.add_argument("--log",
                            help="Optional log file (default: stdout)", default=None)
        args = parser.parse_args()

        reformat_fasta_headers(
            args.accessions,
            args.taxid_mapping,
            args.database_dir,
            args.output_fasta,
            args.output_report,
            log=args.log)
