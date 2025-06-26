import sys
import pandas as pd
import logging
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

# assembly levels in order of quality (best first)
ASSEMBLY_LEVELS = [
    "complete genome",
    "chromosome",
    "scaffold",
    "contig"
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


def get_taxonomic_ranks(cutoff_rank):
    """Return all taxonomic ranks at or lower than the cutoff rank."""
    try:
        idx = TAXONOMIC_RANKS.index(cutoff_rank.lower())
        return TAXONOMIC_RANKS[idx:]
    except ValueError:
        raise ValueError(
            f"'{cutoff_rank}' is not a valid NCBI taxonomic rank.")


def filter_taxonomic_levels(df, max_tax_level):
    """Filter DataFrame to include only rows with taxonomic levels at or below max_tax_level."""
    if max_tax_level:
        max_tax_level = max_tax_level.lower()
        if max_tax_level not in TAXONOMIC_RANKS:
            raise ValueError(
                f"'{max_tax_level}' is not a valid NCBI taxonomic rank.")
        allowed_levels = get_taxonomic_ranks(max_tax_level)
        df['tax_level'] = df['tax_level'].astype(str).str.lower()
        return df[df['tax_level'].isin(allowed_levels)]
    return df


def filter_assembly_levels(df, min_rank):
    """Filter DataFrame to include only rows with assembly levels at or above min_asm_level."""
    df = df[df['asm_level_rank'].notnull()]  # remove invalid assembly levels
    # filter those above the minimum rank
    df = df[df['asm_level_rank'] <= min_rank]
    return df


def keep_top_accessions_per_taxon(df, max_per_taxon):
    """Keep only the top max_per_taxon accessions for each taxon."""
    if max_per_taxon:
        df = df.sort_values(
            by=['taxid', 'asm_level_rank'], ascending=[True, True])
        return df.groupby('taxid').head(max_per_taxon)
    return df


def filter_accessions(input_accessions, output_filtered,
                      max_tax_level='species', max_per_taxon=100,
                      min_asm_level='contig', log=None):

    setup_logging(logfile=log)
    logging.info("Starting accession filtering process.")
    logging.info(f"Input accessions file: {input_accessions}")
    logging.info(f"Output filtered accessions file: {output_filtered}")
    logging.info(f"Maximum taxonomic level: {max_tax_level}")
    logging.info(f"Maximum accessions per taxon: {max_per_taxon}")
    logging.info(f"Minimum assembly level: {min_asm_level}")
    # Load the input accessions
    df = pd.read_csv(input_accessions, sep="\t")
    logging.info("Input accessions loaded successfully.")
    # Filter by taxonomic level
    if max_tax_level:
        size_before = len(df)
        df = filter_taxonomic_levels(df, max_tax_level)
        size_after = len(df)
        logging.info(
            f"Filtered accessions by taxonomic level: {size_before} -> {size_after} entries")

    # Define assembly level order (best first)
    asm_rank = {level: i for i, level in enumerate(ASSEMBLY_LEVELS)}
    df['asm_level_rank'] = df['assembly_level'].str.lower().map(asm_rank)

    # Filter by assembly level
    if min_asm_level:
        size_before = len(df)
        min_rank = asm_rank.get(min_asm_level.lower())
        if min_rank is None:
            raise ValueError(
                f"'{min_asm_level}' is not a valid assembly level.")
        df = filter_assembly_levels(df, min_rank)
        size_after = len(df)
        logging.info(
            f"Filtered accessions by assembly level: {size_before} -> {size_after} entries")

    # Get top accessions per taxid
    if max_per_taxon:
        size_before = len(df)
        df = df.sort_values(
            by=['taxid', 'asm_level_rank'], ascending=[True, True])
        df = df.groupby('taxid').head(max_per_taxon)
        size_after = len(df)
        logging.info(
            f"Keep maximum number per taxon: {size_before} -> {size_after} entries")
    # Save the filtered accessions
    df['accession'].to_csv(output_filtered, sep="\t",
                           index=False, header=False)

    logging.info(f"Filtered accessions saved to {output_filtered}")


if __name__ == "__main__":

    try:
        snakemake
        filter_accessions(
            snakemake.input.accessions,
            snakemake.output.filteredaccessions,
            snakemake.params.max_tax_level,
            snakemake.params.max_per_taxon,
            snakemake.params.min_asm_level,
            log=snakemake.log if snakemake.log else None
        )
    except NameError:

        parser = argparse.ArgumentParser(
            description="Filter NCBI datasets genome accessions."
        )

        parser.add_argument(
            "input_accessions",
            help="Input TSV file with taxid, accession, assembly level, and taxonomic level")
        parser.add_argument(
            "output_filtered", help="Output filtered TSV file")
        parser.add_argument(
            "--log", help="Optional log file (default: stdout)", default=None)
        parser.add_argument(
            "--max_tax_level", type=str, choices=TAXONOMIC_RANKS, default=None,
            help="Maximum taxonomic level to include (default: None, includes all levels)")
        parser.add_argument(
            "--max_per_taxon", type=int, default=None,
            help="Maximum number of accessions per taxon (default: None, includes all accessions)")
        parser.add_argument(
            "--min_asm_level", type=str, choices=ASSEMBLY_LEVELS,
            default=None,
            help="Minimum assembly level to include (default: None, includes all levels)")
        args = parser.parse_args()

        filter_accessions(args.input_accessions, args.output_filtered, args.max_tax_level,
                          args.max_per_taxon, args.min_asm_level, log=args.log)
