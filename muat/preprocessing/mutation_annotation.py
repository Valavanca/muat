# _dmm_annotate_mutations_with_gc_content.py is the original script from the MuAT package.

from pathlib import Path

from collections import deque
import fire
from loguru import logger
from tqdm import tqdm
import pandas as pd
import pybedtools

from util import read_reference, openz


DEFAULT_WINDOW_SIZE = 1001
DEFAULT_BATCH_SIZE = 100000
DEFAULT_LABEL = 'gc'

class MutationAnnotator:
    def __init__(self, window: int, batch_size: int,
                 reference: str, verbose: bool = False):
        self.window = window
        self.batch_size = batch_size
        self.reference = reference
        self.verbose = verbose
        self.ref = read_reference(reference, verbose)

    def calculate_gc_ratio(self, chrom: str, pos: int) -> float:
        """Calculate GC ratio for a given chromosome and position."""
        cpos = max(0, pos - self.window // 2)
        mpos = min(len(self.ref[chrom]) - 1, cpos + self.window)
        buf = deque(self.ref[chrom][cpos:mpos])
        gc = sum(1 for c in buf if c in 'CG')
        at = sum(1 for c in buf if c in 'AT')
        return 1.0 * gc / (gc + at) if gc + at > 0 else 0

    def annotate(self, input_file: str, output_file: str, label: str = DEFAULT_LABEL):
        """Annotate mutations with GC content."""
        with openz(output_file, 'wt') as o, openz(input_file) as f, tqdm(total=self._get_file_line_count(input_file)) as pbar:
            hdr = f.readline().strip().decode("utf-8").split('\t')
            if hdr[0] == 'chrom':
                o.write('{}\t{}\n'.format('\t'.join(hdr), label))

            cchrom = None
            for s in f:
                v = s.strip().decode("utf-8").split('\t')
                chrom, pos = v[0], int(v[1])
                if cchrom != chrom:
                    cchrom = chrom
                gc_ratio = self.calculate_gc_ratio(chrom, pos)
                o.write('{}\t{}\n'.format(s.strip().decode("utf-8"), gc_ratio))
                pbar.update(1)

    def _get_file_line_count(self, file_path: str) -> int:
        """Get the number of lines in a file."""
        with openz(file_path) as f:
            return sum(1 for _ in f)

def annotate_mutations_with_gc_content(
        input_file,
        output_file,
        reference,
        window=DEFAULT_WINDOW_SIZE,
        batch_size=DEFAULT_BATCH_SIZE,
        verbose=False,
        label='gc'
    ):
    """
    Function to annotate mutations with guanine-cytosine content(GC) content. \
    GC-content is the percentage of nitrogenous bases in a DNA or RNA molecule \
    that are either guanine (G) or cytosine (C).
    Original code is `preprocessing/dmm/annotate_mutations_with_gc_content.py`.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        reference (str): Path to the reference genome file.
        window (int): Size of the window for GC content calculation.
        batch_size (int): Batch size for processing.
        verbose (bool): Flag to enable verbose logging.
        label (str): Label for the annotation.

    Returns:
        None
    """
    if verbose:
        logger.info("Starting mutation annotation process")
    annotator = MutationAnnotator(window, batch_size, reference, verbose)
    annotator.annotate(input_file, output_file, label)
    if verbose:
        logger.info("Mutation annotation completed")


def annotate_mutations_with_bed(mutations_tsv, annotations_bed, output_tsv, label='genic'):
    """Annotate a mutation file with the 5th column of a BED file.
    Original script is `preprocessing/dmm/annotate_mutations_with_bed.sh`.

    Args:
        mutations_tsv (str): Path to the mutations TSV file.
        annotations_bed (str): Path to the annotations BED file.
        output_tsv (str): Path to the output TSV file.
        label (str): Column name for the annotation in the output file.

    Returns:
        None
    """
    # Read mutation file
    mutations_df = pd.read_csv(mutations_tsv, sep='\t', low_memory=False)

    # Create a BedTool object from the mutations dataframe
    mutations_bed = pybedtools.BedTool.from_dataframe(mutations_df)

    # Read annotations BED file
    annotations_bed = pybedtools.BedTool(annotations_bed)

    # Perform bedmap operation
    annotated = mutations_bed.map(annotations_bed, c=5, o='mean', null='nan').to_dataframe()

    # Rename the last column and drop extra columns
    annotated.columns = list(mutations_df.columns) + [label]
    annotated = annotated.iloc[:, :-1]

    # Write to output file
    annotated.to_csv(output_tsv, sep='\t', index=False)

    logger.info(f"Annotation completed. Output written to {output_tsv}")


def annotate_with_coding_strand(input_file, output_file, annotation_file, reference_file):
    """
    Annotate a mutation TSV with coding strand information.

    Args:
        input_file (str): Path to the input mutations TSV file.
        output_file (str): Path to the output TSV file.
        annotation_file (str): Path to the BED file with transcript directionality.
        reference_file (str): Path to the reference genome file.
    """
    logger.info("Reading reference genome")
    reference = read_reference(reference_file, verbose=True)

    logger.info("Processing mutations")
    mutations_df = pd.read_csv(input_file, sep='\t', low_memory=False)
    mutations_bed = pybedtools.BedTool.from_dataframe(mutations_df)
    annotations_bed = pybedtools.BedTool(annotation_file)

    annotated = mutations_bed.map(annotations_bed, c=5, o='collapse', null='nan').to_dataframe()
    annotated['strand'] = annotated.apply(lambda row: determine_strand(row, reference), axis=1)

    logger.info("Writing output")
    annotated.to_csv(output_file, sep='\t', index=False)

    logger.info(f"Annotation completed. Output written to {output_file}")

def determine_strand(row, reference):
    """Determine the strand based on reference and overlapping annotations.

    Args:
        row (pandas.Series): A row from the dataframe.
        reference (dict): A dictionary containing the reference genome sequences.

    Returns:
        str: The strand information.
    """
    chrom, pos = row['chrom'], row['start']
    base = reference[chrom][pos]
    strands = row[4]  # Assuming the strand information is in the 5th column

    if '+' in strands and '-' in strands:
        return '?'
    elif '+' in strands:
        return '+' if base in ['C', 'T'] else '-'
    elif '-' in strands:
        return '-' if base in ['C', 'T'] else '+'
    else:
        return '='


if __name__ == '__main__':
    fire.Fire(
        {
            'annotate_mutations_with_gc_content': annotate_mutations_with_gc_content,
            'annotate_mutations_with_bed': annotate_mutations_with_bed,
            'annotate_with_coding_strand': annotate_with_coding_strand
        }
    )