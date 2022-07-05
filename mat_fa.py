from Bio import SeqIO
from pathlib import Path
import gzip
from desc_fa import find_in_description


def parse_file(filepath, ext, gz, func_name, **kwargs):
    """
    Accepts a file with nucleotide sequences and parses it.
    Then applies a specified function to each record and writes the modified records to a file.
    Returns the path to the created file.
    :param filepath: path to the file with nucleotide sequences.
    :param ext: extension of the provided file: fasta or fastq. The output file will have the same extension,
    except for the func_name == translate: in this case the output format is fasta in order to prevent
    possible letter_annotations conflicts.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :param func_name: name of the function to apply to each record in the file.
    :param kwargs: arguments of the function applied if present.
    :return: path to the created file.
    """
    out_dir = Path(filepath).parent
    out_stem = Path(filepath).stem.split('.')[0]
    if func_name == translate:
        out_file = out_stem + '_' + func_name.__name__ + '.fasta'
    else:
        out_file = out_stem + '_' + func_name.__name__ + '.' + ext
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt", encoding='utf-8') as out_handle:
            for record in records:
                func_name(record, **kwargs)
                if func_name == translate:
                    SeqIO.write(record, out_handle, 'fasta')
                else:
                    SeqIO.write(record, out_handle, ext)
    else:
        with open(out_path, "w", encoding='utf-8') as out_handle:
            records = SeqIO.parse(filepath, ext)
            for record in records:
                func_name(record, **kwargs)
                if func_name == translate:
                    SeqIO.write(record, out_handle, 'fasta')
                else:
                    SeqIO.write(record, out_handle, ext)
    return out_path


def complement(record):
    record.seq = record.seq.complement()


def reverse_complement(record):
    record.seq = record.seq.reverse_complement()


def transcribe(record):
    record.seq = record.seq.transcribe()


def back_transcribe(record):
    record.seq = record.seq.back_transcribe()


def translate(record, table=1, stop_symbol='*', to_stop=False, cds=False, gap=None):
    """
    Performs translation and provides a protein sequence for a nucleotide sequence.
    :param record: a record containing a nucleotide sequence to translate.
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
    record.letter_annotations = {}
    while len(record.seq) % 3 != 0:
        record.seq = record.seq[:-1]
    record.seq = record.seq.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap)


if __name__ == "__main__":

    # Example 1. Find all DNA sequences related to chromosome 7, then transcribe and translate them.
    chr7 = find_in_description('./test_files/Sample_1/Sample_1.fasta', 'chromosome 7', 'fasta', gz=False)
    ch7_transcribed = parse_file(chr7, 'fasta', False, transcribe)
    parse_file(ch7_transcribed, 'fasta', False, translate)

    # Example 2. Make reverse complements for DNA sequences, transcribe and translate them.
    rev_com = parse_file('./test_files/Sample_2/Sample_2.fastq.gz', 'fastq', True, reverse_complement)
    rc_transcribed = parse_file(rev_com, 'fastq', True, transcribe)
    parse_file(rc_transcribed, 'fastq', True, translate, table=2)

    # Example 3. Back transcribe RNA sequences and translate DNA.
    tr = parse_file('./test_files/Sample_3/Sample_3_RNA.fastq.gz', 'fastq', True, back_transcribe)
    parse_file(tr, 'fastq', True, translate)
