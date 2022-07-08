from Bio import SeqIO
from pathlib import Path
import gzip


def _get_out_path(filepath, out_str, ext, gz):
    out_dir = Path(filepath).parent
    if gz:
        out_file = out_str.replace(" ", "_") + "." + ext + ".gz"
    else:
        out_file = out_str.replace(" ", "_") + "." + ext
    return Path(out_dir, out_file)


def find_in_description(filepath, str_to_find):
    ext_dot = Path(filepath).suffix
    ext = ext_dot[1:]
    out_str = Path(filepath).stem.split('.')[0] + "_" + str_to_find
    out_path = _get_out_path(filepath=filepath, out_str=out_str, ext=ext, gz=False)
    matched_recs = (rec for rec in SeqIO.parse(filepath, ext) if str_to_find in rec.description)
    SeqIO.write(matched_recs, out_path, ext)
    return out_path


def find_id_description_gz(filepath, str_to_find):
    ext_dot = Path(filepath).suffixes
    ext = ext_dot[-2][1:]
    out_str = Path(filepath).stem.split('.')[0] + "_" + str_to_find
    out_path = _get_out_path(filepath=filepath, out_str=out_str, ext=ext, gz=True)
    matched_recs = (rec for rec in SeqIO.parse(gzip.open(filepath, mode='rt', encoding='utf-8'), ext) if str_to_find in rec.description)
    with gzip.open(out_path, "wt") as out_handle:
        SeqIO.write(matched_recs, out_handle, ext)
    return out_path
