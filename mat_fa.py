from Bio import SeqIO
from pathlib import Path
import gzip
from search_fa import find_in_description
from sort_fa import sort_by_length
from qual_fastq import high_quality_gz


def _get_out_path(filepath, func_name, ext, gz):
    out_dir = Path(filepath).parent
    out_stem = Path(filepath).stem.split('.')[0]
    if gz:
        if func_name == translate:
            out_file = f"{out_stem}_{func_name.__name__}.fasta.gz"
        else:
            out_file = f"{out_stem}_{func_name.__name__}.{ext}.gz"
    else:
        if func_name == translate:
            out_file = f"{out_stem}_{func_name.__name__}.fasta"
        else:
            out_file = f"{out_stem}_{func_name.__name__}.{ext}"
    return Path(out_dir, out_file)


def parse_file(filepath, func_name, **kwargs):
    """
    Accepts a fasta/fastq file with nucleotide sequences, applies the specified function
    to each record and writes the modified record to a file. The function applied represents
    a matrix process, e.g.: the function transcribe() performs transcription.
    The output file has the same format as the input one, except for the func_name == translate:
    in this case the output file is fasta. Returns a path to the created file.
    :param filepath: path to the input file.
    :param func_name: the function to apply.
    :param kwargs: the function's arguments if any.
    :return: path to the output file.
    """
    ext_dot = Path(filepath).suffix
    ext = ext_dot[1:]
    out_path = _get_out_path(filepath=filepath, func_name=func_name, ext=ext, gz=False)
    records = func_name(SeqIO.parse(filepath, ext), **kwargs)
    if func_name == translate:
        SeqIO.write(records, out_path, "fasta")
    else:
        SeqIO.write(records, out_path, ext)
    return out_path


def parse_gz_file(filepath, func_name, **kwargs):
    """
    Accepts a gzipped fasta/fastq file with nucleotide sequences, applies the specified function
    to each record and writes the modified record to a gzipped file. The function applied represents
    a matrix process, e.g.: the function transcribe() performs transcription.
    The output file has the same format as the input one, except for the func_name == translate:
    in this case the output file is fasta. Returns a path to the created file.
    :param filepath: path to the input file.
    :param func_name: the function to apply.
    :param kwargs: the function's arguments if any.
    :return: path to the output file.
    """
    ext_dot = Path(filepath).suffixes
    ext = ext_dot[-2][1:]
    out_path = _get_out_path(filepath=filepath, func_name=func_name, ext=ext, gz=True)
    records = func_name(SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext), **kwargs)
    if func_name == translate:
        with gzip.open(out_path, "wt") as out_handle:
            SeqIO.write(records, out_handle, "fasta")
    else:
        with gzip.open(out_path, "wt") as out_handle:
            SeqIO.write(records, out_handle, ext)
    return out_path


def complement(records):
    for record in records:
        record.seq = record.seq.complement()
        yield record


def reverse_complement(records):
    for record in records:
        record.seq = record.seq.reverse_complement()
        yield record


def transcribe(records):
    for record in records:
        record.seq = record.seq.transcribe()
        yield record


def back_transcribe(records):
    for record in records:
        record.seq = record.seq.back_transcribe()
        yield record


def translate(records, table=1, stop_symbol='*', to_stop=False, cds=False, gap=None):
    """
    Performs translation and provides a protein sequence for a nucleotide sequence.
    :param records: records with nucleotide sequences to translate.
    :param table: translation table, string or integer, according to
    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. The Standard Code by default.
    :param stop_symbol: single character string, what to use for terminators. This defaults to the asterisk, “*”.
    :param to_stop: boolean, True if you want to stop translation after the first stop codon,
    False otherwise. Default is False.
    :param cds: Boolean, indicates this is a complete CDS.
    If True, this checks the sequence starts with a valid alternative start codon
    (which will be translated as methionine, M), that the sequence length is
    a multiple of three, and that there is a single in frame stop codon at the end
    (this will be excluded from the protein sequence, regardless of the to_stop option).
    If these tests fail, an exception is raised.
    :param gap: Single character string to denote symbol used for gaps.
    It will try to guess the gap character from the alphabet.
    """
    for record in records:
        record.letter_annotations = {}
        while len(record.seq) % 3 != 0:
            record.seq = record.seq[:-1]
        record.seq = record.seq.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap)
        yield record


if __name__ == "__main__":

    # Example 1. Find all DNA sequences related to chromosome 7, then transcribe and translate them.
    chr7 = find_in_description('./test_files/Sample_1/Sample_1.fasta', 'chromosome 7')
    ch7_transcribed = parse_file(chr7, transcribe)
    parse_file(ch7_transcribed, translate)

    # Example 2. Make reverse complements for DNA sequences, transcribe and translate them.
    qual_20 = high_quality_gz('./test_files/Sample_2/Sample_2.fastq.gz', 20)
    rev_com = parse_gz_file(qual_20, reverse_complement)
    rc_transcribed = parse_gz_file(rev_com, transcribe)
    parse_gz_file(rc_transcribed, translate, table=2)

    # Example 3. Sorting by lengths.
    sort_by_length('./test_files/Sample_2/Sample_2.fastq')

    # Example 4. Back transcribe RNA sequences and translate DNA.
    back_tr = parse_gz_file('./test_files/Sample_3/Sample_3_RNA.fastq.gz', back_transcribe)
    parse_gz_file(back_tr, translate)

