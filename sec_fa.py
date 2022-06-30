from Bio import SeqIO
from pathlib import Path
import gzip
from desc_fa import find_in_description


def _get_dir_stem(filepath):
    in_file = Path(filepath)
    return [in_file.parent, in_file.stem.split('.')[0]]


def complement(filepath, ext, gz, reverse=False):
    """
    Writes a complement or a reverse complement for all nucleotide sequences from a given file.
    :param filepath: path to the file with sequences to write complements for.
    :param ext: xtension of the provided file: fasta or fastq. The output file will have the same extension.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :param reverse: boolean, if True writes reverse complements for the provided sequences.
    :return:
    """
    out_dir, out_stem = _get_dir_stem(filepath)
    if reverse:
        out_file = out_stem + '_reverse_complement.' + ext
    else:
        out_file = out_stem + '_complement.' + ext
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt") as out_handle:
            for record in records:
                if reverse:
                    record.seq = record.seq.reverse_complement()
                else:
                    record.seq = record.seq.complement()
                SeqIO.write(record, out_handle, ext)
    else:
        with open(out_path, "w") as out_handle:
            records = SeqIO.parse(filepath, ext)
            for record in records:
                if reverse:
                    record.seq = record.seq.reverse_complement()
                else:
                    record.seq = record.seq.complement()
                SeqIO.write(record, out_handle, ext)
    return out_path


def transcribe_dna(filepath, ext, gz):
    """
    Performs transcription of DNA sequence(s) from the provided file.
    Writes the records with transcribed sequences to a file.
    :param filepath: path to the file with DNA sequence(s) to transcribe.
    :param ext: extension of the provided file: fasta or fastq. The output file will have the same extension.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :return: path to the created file.
    """
    out_dir, out_stem = _get_dir_stem(filepath)
    out_file = out_stem + '_transcribed.' + ext
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt") as out_handle:
            for record in records:
                record.seq = record.seq.transcribe()
                SeqIO.write(record, out_handle, ext)
    else:
        with open(out_path, "w") as out_handle:
            records = SeqIO.parse(filepath, ext)
            for record in records:
                record.seq = record.seq.transcribe()
                SeqIO.write(record, out_handle, ext)
    return out_path


def back_transcribe_rna(filepath, ext, gz):
    """
    Performs back transcription of RNA sequence(s) from the provided file.
    Writes the records with back transcribed sequences to a file.
    :param filepath: path to the file with RNA sequence(s) to transcribe back.
    :param ext: extension of the provided file: fasta or fastq. The output file will have the same extension.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :return: path to the created file.
    """
    out_dir, out_stem = _get_dir_stem(filepath)
    out_file = out_stem + '_back_transcribed.' + ext
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt") as out_handle:
            for record in records:
                record.seq = record.seq.back_transcribe()
                SeqIO.write(record, out_handle, ext)
    else:
        with open(out_path, "w") as out_handle:
            records = SeqIO.parse(filepath, ext)
            for record in records:
                record.seq = record.seq.back_transcribe()
                SeqIO.write(record, out_handle, ext)
    return out_path


def translate_dna_rna(filepath, ext, gz, to_stop=False, table=1):
    """
    Performs translation of RNA (or DNA) sequence(s) from the provided fil.
    Writes the record with translated sequences to a fasta file.
    :param filepath: path to the file with RNA (DNA) sequence(s) to translate.
    :param ext: extension of the provided file: fasta or fastq. The output file will be fasta in oreder
    to prevent possible letter_annotations conflicts.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :param to_stop: boolean, Ture if you want to stop translation after the first stop codon, False otherwise. Default is False.
    :param table: translation table, string or integer, according to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
    The Standard Code by default
    :return: path to the created file.
    """
    out_dir, out_stem = _get_dir_stem(filepath)
    out_file = out_stem + '_translated.fasta'
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt") as out_handle:
            for record in records:
                record.letter_annotations = {}
                while len(record.seq) % 3 != 0:
                    record.seq = record.seq[:-1]
                record.seq = record.seq.translate(to_stop=to_stop, table=table)
                SeqIO.write(record, out_handle, 'fasta')
    else:
        with open(out_path, "w") as out_handle:
            records = SeqIO.parse(filepath, ext)
            for record in records:
                record.letter_annotations = {}
                while len(record.seq) % 3 != 0:
                    record.seq = record.seq[:-1]
                record.seq = record.seq.translate(to_stop=to_stop, table=table)
                SeqIO.write(record, out_handle, 'fasta')
    return out_path


if __name__ == "__main__":

    # Example 1. Find all DNA sequences related to chromosome 4, then transcribe and translate the.
    chr4 = find_in_description('./Sample_1/Sample_1.fasta', 'chromosome 4', 'fasta', gz=False)
    ch4_transcribed = transcribe_dna(chr4, 'fasta', gz=False)
    translate_dna_rna(ch4_transcribed, 'fasta', gz=False)

    # Example 2. Make reverse complements, then transcribe and then translate the complemen sequences.
    rev_comp = complement('./Sample_2/Sample_2.fastq.gz', 'fastq', gz=True, reverse=True)
    rc_transcribed = transcribe_dna(rev_comp, 'fastq', gz=True)
    translate_dna_rna(rc_transcribed, 'fastq', gz=True)

    # Example 3. Back transcription for RNA sequences.
    back_transcribe_rna('./Sample_3/Sample_3_RNA.fastq.gz', 'fastq', gz=True)





