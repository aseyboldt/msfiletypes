import xml.etree.ElementTree as ET
import tables as tb
from scipy import sparse
import argparse
import contextlib


class Proteins(tb.IsDescription):
    id = tb.Int64Col(pos=0)
    accession = tb.StringCol(200, pos=1)
    score = tb.Float64Col(pos=2)


class PeptideHits(tb.IsDescription):
    pos_id = tb.Int64Col(pos=0)
    mz = tb.Float64Col(pos=1)
    rt = tb.Float64Col(pos=2)
    charge = tb.Int64Col(pos=3)
    score = tb.Float64Col(pos=4)
    sequence = tb.StringCol(500, pos=5)
    aa_before = tb.StringCol(1, pos=6)
    aa_after = tb.StringCol(1, pos=7)
    e_value = tb.Float64Col(pos=8)


def write_run(search_params, ident_run, hdf_file, where="/", filters=None):
    """ Store the data from one identification run in a hdf5 group.

    Parameters
    ----------
    search_params: ElementTree.Element
        A 'SearchParameters' section of idXML
    ident_run: ElementTree.Element
        The corresponding 'IdentificationRun' section.
    hdf_file: tables.File
        Opened HDF5 file from pyTables
    where: str or group
        The base element as in tables.create_table
    filters: tables.Filter
        Filters to use in HDF5. See tables.Filter
    """
    if not search_params.tag == 'SearchParameters':
        raise ValueError('Not a SearchParameters tag: %s' % search_params)
    if not ident_run.tag == 'IdentificationRun':
        raise ValueError('Not a IdentificationRun tag: %s' % ident_run)

    proteins = hdf_file.create_table(
        where, 'proteins', description=Proteins, filters=filters
    )

    prot_ident_el = ident_run.find('ProteinIdentification')
    assert prot_ident_el is not None
    for i, prot in enumerate(prot_ident_el.findall('ProteinHit')):
        prot_id = i
        accession = prot.get('accession')
        score = float(prot.get('score'))
        proteins.append([(prot_id, accession, score)])

    peptide_hits = hdf_file.create_table(
        where, 'peptide_hits', filters=filters, description=PeptideHits
    )

    hit_to_prot = set()

    peptide_hit_num = 0
    peptide_hit_groups = ident_run.findall('PeptideIdentification')
    for i, peptide_hit_group in enumerate(peptide_hit_groups):
        pos_id = i
        mz = float(peptide_hit_group.get('MZ'))
        rt = float(peptide_hit_group.get('RT'))
        for peptide_hit in peptide_hit_group.findall('PeptideHit'):
            seq = peptide_hit.get('sequence')
            charge = int(peptide_hit.get('charge'))
            score = float(peptide_hit.get('score'))
            aa_before = peptide_hit.get('aa_before')
            aa_after = peptide_hit.get('aa_after')
            e_value = peptide_hit[0].get('value')
            prot_refs = peptide_hit.get('protein_refs').split(' ')
            for ref in prot_refs:
                hit_to_prot.add((peptide_hit_num, int(ref[3:])))
            peptide_hits.append([
                (pos_id, mz, rt, charge, score, seq, aa_before,
                 aa_after, e_value)
            ])
            peptide_hit_num += 1

    m = sparse.dok_matrix((len(peptide_hits), len(proteins)), dtype=bool)
    for i, j in hit_to_prot:
        m[i, j] = True

    m = m.tocsr()
    hdf_file.create_array(where, 'identification_indptr', obj=m.indptr)
    hdf_file.create_array(where, 'identification_indices', obj=m.indices)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an OpenMS idXML file to a hdf5 format"
    )

    parser.add_argument('input', help='input idXML file', type=str)
    parser.add_argument('output', help='Output path')
    parser.add_argument('--compress', choices=['blosc', 'zlib', 'lzo'],
                        default='zlib')
    parser.add_argument('--compression-level', choices=list(range(1, 10)),
                        default=1)
    return parser.parse_args()


def main():
    args = parse_args()
    infile = ET.parse(args.input)
    root = infile.getroot()

    filters = tb.Filters(complevel=args.compression_level,
                         complib=args.compress, fletcher32=True)
    f = tb.open_file(args.output, 'w')
    with contextlib.closing(f):
        run_group = f.create_group('/', 'run000')
        write_run(root[0], root[1], f, where=run_group, filters=filters)

if __name__ == '__main__':
    main()
