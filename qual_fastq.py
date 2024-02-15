from Bio import SeqIO
from pathlib import Path
import gzip


def _get_out_path(filepath, out_str, ext, gz):
    out_dir = Path(filepath).parent
    out_stem = Path(filepath).stem.split('.')[0]
    if gz:
        out_file = f"{out_stem}_{out_str}.{ext}.gz"
    else:
        out_file = f"{out_stem}_{out_str}.{ext}"
    return Path(out_dir, out_file)


def high_quality(filepath, qual):
    """
    Accepts a fastq file and searches for reads having each symbol of the sequence
    with the quality not less than specified.
    :param filepath: path to the fastq file.
    :param qual: int, the quality to filter the reads by. All reads in the output file will have
    this quality or higher.
    :return: the path to the output file.
    """
    ext_dot = Path(filepath).suffix
    ext = ext_dot[1:]
    out_str = "min_qual_" + str(qual)
    out_path = _get_out_path(filepath=filepath, out_str=out_str, ext=ext, gz=False)
    matched_recs = (
        rec for rec in SeqIO.parse(filepath, ext)
        if min(rec.letter_annotations["phred_quality"]) >= qual
    )
    SeqIO.write(matched_recs, out_path, ext)
    return out_path


def high_quality_gz(filepath, qual):
    """
    Accepts a gzipped fastq file and searches for reads having each symbol of the sequence
    with the quality not less than specified.
    :param filepath: path to the fastq file.
    :param qual: int, the quality to filter the reads by. All reads in the output file will have
    this quality or higher.
    :return: the path to the output file.
    """
    ext_dot = Path(filepath).suffixes
    ext = ext_dot[-2][1:]
    out_str = "min_qual_" + str(qual)
    out_path = _get_out_path(filepath=filepath, out_str=out_str, ext=ext, gz=True)
    matched_recs = (
        rec for rec in SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext)
        if min(rec.letter_annotations["phred_quality"]) >= qual
    )
    with gzip.open(out_path, "wt") as out_handle:
        SeqIO.write(matched_recs, out_handle, ext)
    return out_path
