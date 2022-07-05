from Bio import SeqIO
from pathlib import Path
import gzip


def find_in_description(filepath, str_to_find, ext, gz):
    """
    Searches for a specified string in record.description in all records of the provided file.
    Writes all matched records to a file.
    :param filepath: path to the file with records to analyse.
    :param str_to_find: string to fins in the record description.
    :param ext: extension of the provided file. The output file will have the same extension.
    Though initially this function was intended for fast and fastq formats, you can choose any
    file format that has record.description after implementation of Biopython function SeqIO.parse.
    See more at https://biopython.org/wiki/SeqIO. Keep in mind that for formats other than fasta of fastq
    the function's behaviour can be unpredicted.
    :param gz: boolean, True if the file is gzipped, False otherwise. The output file will have the same gz value.
    :return: path to the created file.
    """
    out_dir = Path(filepath).parent
    out_file = Path(filepath).stem.split('.')[0] + "_" + str_to_find.replace(" ", "_") + "." + ext
    out_path = Path(out_dir, out_file)
    if gz:
        out_path = Path(out_dir, out_file + ".gz")
        records = SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        with gzip.open(out_path, "wt", encoding='utf-8') as out_handle:
            for record in records:
                if str_to_find in record.description:
                    SeqIO.write(record, out_handle, ext)
    else:
        records = SeqIO.parse(filepath, ext)
        with open(out_path, "w", encoding='utf-8') as out_handle:
            for record in records:
                if str_to_find in record.description:
                    SeqIO.write(record, out_handle, ext)
    return out_path
