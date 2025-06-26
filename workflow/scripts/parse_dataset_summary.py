import sys
import json
import pandas as pd
import numpy as np
from collections import Counter
import logging
from ete3 import NCBITaxa
import argparse

# Ordered list of standard NCBI taxonomic ranks (top to bottom)
TAXONOMIC_RANKS = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain"
]


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


def read_json(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)


def parse_json(data, rollup_tax_level=None):
    total_assemblies = data.get('total_count', 0)
    reports = data.get('reports', [])
    rows = []
    for assembly in reports:
        taxid = assembly.get('organism', {}).get('tax_id')
        accession = assembly.get('accession')
        asm_level = assembly.get('assembly_info', {}).get('assembly_level')
        rows.append([taxid, accession, asm_level])
    df = pd.DataFrame(rows, columns=['taxid', 'accession', 'assembly_level'])
    return {'data': df, 'total_assemblies': total_assemblies}


def get_rollup_taxonomic_level(df, tax_level, ncbi):
    """Return the taxid at provided level in the lineage of the given taxid."""
    try:
        idx = TAXONOMIC_RANKS.index(tax_level.lower())
    except ValueError:
        raise ValueError(f"'{tax_level}' is not a valid NCBI taxonomic rank.")

    ranks = ncbi.get_rank(df['taxid'].unique().tolist())
    df['tax_level'] = df['taxid'].map(ranks)
    roll_up_taxa = df[df['tax_level'] != tax_level.lower()]['taxid'].unique()
    roll_up_map = {}
    for taxid in roll_up_taxa:
        lineage = ncbi.get_lineage(taxid)
        ranks = {v: k for k, v in ncbi.get_rank(lineage).items()}
        for rank in TAXONOMIC_RANKS[idx:0:-1]:
            if rank in ranks.keys():
                roll_up_map[taxid] = np.int64(ranks[rank])
                break
    df['taxid'] = df['taxid'].apply(lambda x: roll_up_map.get(x, x))
    return df.drop('tax_level', axis=1, errors='ignore')


def parse_dataset_summary(json_file, accession_list, summary, min_tax_level=None, log=None):
    setup_logging(logfile=log)

    logging.info("Starting to parse dataset summary JSON file.")
    logging.info(f"Input JSON file: {json_file}")
    logging.info(f"Output TSV file: {accession_list}")
    logging.info(f"Report file: {summary}")
    logging.info(
        f"Roll up taxa to minimum taxonomic level: {min_tax_level if min_tax_level else 'None'}")
    try:
        data = read_json(json_file)
        logging.info("JSON file read successfully.")
    except Exception as e:
        logging.error(f"Error reading JSON file: {e}")
        return

    result = parse_json(data)
    logging.info("JSON file parsed successfully.")
    logging.info(f"Total assemblies found: {result['total_assemblies']}")

    ncbi = NCBITaxa()

    df = result['data']
    if min_tax_level:
        df = get_rollup_taxonomic_level(df, min_tax_level, ncbi)
    taxid_list = df['taxid'].unique().tolist()
    logging.info(f"Total unique taxa: {len(taxid_list)}")

    ranks = ncbi.get_rank(taxid_list)

    df['tax_level'] = df['taxid'].map(ranks)
    taxa_names = ncbi.get_taxid_translator(taxid_list)
    df['tax_name'] = df['taxid'].map(taxa_names)
    df.to_csv(accession_list, sep="\t", index=False)
    logging.info(f"Data written to TSV file: {accession_list}")

    # Create summary report
    counts = pd.DataFrame.from_dict(
        df['taxid'].value_counts().to_dict(), orient='index', columns=['count'])
    counts.index.name = 'taxid'
    counts.reset_index(inplace=True)
    counts['tax_level'] = counts['taxid'].map(ranks)
    counts['tax_name'] = counts['taxid'].map(taxa_names)
    counts = counts.sort_values(by='count', ascending=False)
    counts = counts[['taxid', 'tax_name', 'tax_level', 'count']]
    logging.info(
        f"Taxonomic levels in output: {counts['tax_level'].unique()}")
    counts.to_csv(summary, sep="\t", index=False)
    logging.info(f"Summary report written to: {summary}")
    logging.info("Parsing completed successfully.")


if __name__ == "__main__":

    try:
        snakemake
        parse_dataset_summary(
            snakemake.input.json,
            snakemake.output.accessions,
            snakemake.output.counts,
            snakemake.params.min_tax_level,
            log=snakemake.log if snakemake.log else None
        )
    except NameError:

        parser = argparse.ArgumentParser(
            description="Parse NCBI datasets genome summary JSON to TSV."
        )

        parser.add_argument(
            "input_json",
            help="Input JSON file from datasets genome summary")
        parser.add_argument(
            "output_tsv",
            help="Output TSV file with taxid, accession, assembly level, and taxonomic level")
        parser.add_argument(
            "report_file",
            help="Summary report of number of assemblies per taxon")
        parser.add_argument(
            "--min-tax-level",
            choices=TAXONOMIC_RANKS, default=None,
            help="Minimum taxonomic level to include in the report, lower will be rolled up (default: all)")
        parser.add_argument(
            "--log", help="Optional log file (default: stdout)", default=None)

        args = parser.parse_args()

        parse_dataset_summary(args.input_json, args.output_tsv,
                              args.report_file, args.min_tax_level, log=args.log)
