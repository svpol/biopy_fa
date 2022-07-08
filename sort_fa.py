from Bio import SeqIO
from pathlib import Path


def _get_out_path(filepath, out_str, ext):
    out_dir = Path(filepath).parent
    out_stem = Path(filepath).stem.split('.')[0]
    out_file = out_stem + "_" + out_str + "." + ext
    return Path(out_dir, out_file)


def sort_by_id(filepath):
    """
    Accepts a fasta/fastq file and sorts them by record.id.
    Writes the sorted records to the file of the same format.
    Returns the path to the output file.
    :param filepath: path to the input file.
    :return: path to the output file.
    """
    ext_dot = Path(filepath).suffix
    ext = ext_dot[1:]
    out_path = _get_out_path(filepath=filepath, out_str='sorted_by_id', ext=ext)
    ids = sorted(rec.id for rec in SeqIO.parse(filepath, ext))
    record_index = SeqIO.index(filepath, ext)
    records = (record_index[id] for id in ids)
    SeqIO.write(records, out_path, ext)
    return out_path


def sort_by_length(filepath):
    """
    Accepts a fasta/fastq file and sorts them by sequence length.
    Writes the sorted records to the file of the same format.
    Returns the path to the output file.
    :param filepath: path to the input file.
    :return: path to the output file.
    """
    ext_dot = Path(filepath).suffix
    ext = str(ext_dot)[1:]
    out_path = _get_out_path(filepath=filepath, out_str='sorted_by_len', ext=ext)
    len_ids = sorted((len(rec), rec.id) for rec in SeqIO.parse(filepath, ext))
    ids = reversed([id for (length, id) in len_ids])
    del len_ids
    record_index = SeqIO.index(filepath, ext)
    records = (record_index[id] for id in ids)
    SeqIO.write(records, out_path, ext)
    return out_path
