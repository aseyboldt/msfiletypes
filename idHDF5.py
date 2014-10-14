import xml.etree.ElementTree as ET
import tables as tb
import argparse
import contextlib
import numpy as np


class Proteins(tb.IsDescription):
    """ Corresponds to a ProteinHit in idXML. Does not store UserParams """
    id = tb.Int64Col(pos=0)
    accession = tb.StringCol(200, pos=1)
    score = tb.Float64Col(pos=2)


class PeptideHits(tb.IsDescription):
    """ Corresponds to PeptideHit in idXML.

    Peptide hits that share the same PeptideIdentification section (they
    could explain the same spectrum) share their pos_id"""
    pos_id = tb.Int64Col(pos=0)
    mz = tb.Float64Col(pos=1)
    rt = tb.Float64Col(pos=2)
    charge = tb.Int64Col(pos=3)
    score = tb.Float64Col(pos=4)
    sequence = tb.StringCol(500, pos=5)
    aa_before = tb.StringCol(1, pos=6)
    aa_after = tb.StringCol(1, pos=7)
    e_value = tb.Float64Col(pos=8)


def write_proteins(iterparser, start_tag, table):

    def read_protein_hit(elem, buffer):
        buffer['accession'] = elem.get('accession')
        buffer['score'] = float(elem.get('score'))
        buffer['id'] = prot_count

    buffer_size = 1024
    buffer_array = np.empty((buffer_size, ), table.dtype)

    prot_count = 0
    buffer_index = 0
    for event, elem in iterparser:
        if event == 'start':
            continue
        if elem.tag == 'ProteinIdentification':
            table.append(buffer_array[:buffer_index])
            start_tag.clear()
            return
        elif elem.tag == 'ProteinHit':
            read_protein_hit(elem, buffer_array[buffer_index])
            buffer_index += 1
            if buffer_index == buffer_size:
                table.append(buffer_array)
                buffer_index = 0
            prot_count += 1
            start_tag.clear()
        elif elem.tag == 'UserParam':
            pass
        else:
            raise ValueError('Got unexpected xml tag: %s' % elem.tag)


def write_peptides(iterparser, start_elem, peptide_hits, indptr, indices):

    def read_peptide_ident(elem, pos_id):
        assert elem.tag == 'PeptideIdentification'
        mz = float(elem.get('MZ'))
        rt = float(elem.get('RT'))
        protein_refs = []
        rows = []
        for peptide_hit in elem.findall('PeptideHit'):
            seq = peptide_hit.get('sequence')
            charge = int(peptide_hit.get('charge'))
            score = float(peptide_hit.get('score'))
            aa_before = peptide_hit.get('aa_before')
            aa_after = peptide_hit.get('aa_after')
            e_value_elem = [
                el for el in peptide_hit
                if el.get('name') == 'E-Value'
            ]
            if len(e_value_elem) > 0:
                assert len(e_value_elem) == 1
                e_value = float(e_value_elem[0].get('value'))
            else:
                e_value = float('nan')
            prot_refs = peptide_hit.get('protein_refs').split(' ')
            protein_refs.append([int(ref[3:]) for ref in prot_refs])
            rows.append((pos_id, mz, rt, charge, score, seq, aa_before,
                         aa_after, e_value))
        return rows, protein_refs

    pos_id = 0
    for event, elem in iterparser:
        if event == 'start':
            continue
        if elem.tag == 'PeptideIdentification':
            rows, idxs = read_peptide_ident(elem, pos_id)
            start_elem.clear()
            peptide_hits.append(rows)
            for idx in idxs:
                indptr.append([len(indices)])
                indices.append(idx)
            pos_id += 1
            start_elem.clear()
        if elem.tag == 'IdentificationRun':
            indptr.append([len(indices)])
            start_elem.clear()
            return


def write_run(iterparser, start_tag, hdf_file, where="/", filters=None):
    """ Store the data from one identification run in a hdf5 group.

    Parameters
    ----------
    search_params: iterator like ElementTree.iterparse
        This must be in a state after ('start', 'IdentificationRun') has been
        called.
    start_tag: ElementTree.Element
        The start tag element of IdentificationRun
    hdf_file: tables.File
        Opened HDF5 file from pyTables
    where: str or group
        The base element as in tables.create_table
    filters: tables.Filter
        Filters to use in HDF5. See tables.Filter
    """
    if not start_tag.tag == 'IdentificationRun':
        raise ValueError('Not a IdentificationRun tag: %s' % start_tag.tag)

    proteins = hdf_file.create_table(
        where, 'proteins', description=Proteins, filters=filters
    )

    peptide_hits = hdf_file.create_table(
        where, 'peptide_hits', filters=filters, description=PeptideHits
    )

    indptr = hdf_file.create_earray(
        where, 'identification_indptr', tb.Int64Atom(), shape=(0,)
    )
    indices = hdf_file.create_earray(
        where, 'identification_indices', tb.Int64Atom(), shape=(0,)
    )

    for event, elem in iterparser:
        if event == 'start' and elem.tag == 'ProteinIdentification':
            write_proteins(iterparser, elem, proteins)
            elem.clear()
        elif event == 'start' and elem.tag == 'PeptideIdentification':
            write_peptides(iterparser, elem, peptide_hits, indptr, indices)
            elem.clear()
            return
        else:
            raise ValueError('Reached unexpected tag: %s' % elem.tag)


def read_idXML(infile, storage, where='/', no_run_group=False, filters=None):
    """ Read itXML and store the arrays in hdf storage.

    Parameters
    ----------
    infile : path or file
        Input idXML file to read
    storage : path or pytables.File instance
        The HDF5 file to store the tables in. If a path is specified the
        file will be opened in append mode.
    where : str or pytables.Group
        The path in the HDF5 file to use
    no_run_group : bool
        If True, this function will create a seperate group for each run
        in the idXML file.
    filters : pytables.Filter
        pytables.Filter instance (e.g. to specify compression)
    """
    def read_groups(storage_file):
        if where not in storage_file:
            storage_file.create_group('/', where.lstrip('/'))
        run_id = 0
        for event, elem in infile:
            if event == 'start' and elem.tag == 'IdentificationRun':
                if not no_run_group:
                    run_group = f.create_group(
                        where, 'run{:03}'.format(run_id)
                    )
                else:
                    run_group = f.get_node(where)
                write_run(infile, elem, f, where=run_group, filters=filters)
                elem.clear()
                run_id += 1

    infile = ET.iterparse(infile, ['start', 'end'])
    if isinstance(storage, tb.File):
        read_groups(storage)
    else:
        f = tb.open_file(storage, 'a')
        with contextlib.closing(f):
            read_groups(f)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an OpenMS idXML file to a hdf5 format"
    )

    parser.add_argument('input', help='input idXML file', type=str)
    parser.add_argument('output', help='Output path')
    parser.add_argument('--compress', choices=['blosc', 'zlib', 'lzo'],
                        default='zlib')
    parser.add_argument('--compression-level', choices=list(range(1, 10)),
                        default=1, type=int)
    parser.add_argument('--where', help='hdf5 destination inside output',
                        default='/')
    parser.add_argument('--no-run-group', action='store_false', default=True,
                        help='Write the first IdentificationRun only and ' +
                        'do not create a group for the run')
    return parser.parse_args()


def main():
    args = parse_args()
    filters = tb.Filters(complevel=args.compression_level,
                         complib=args.compress, fletcher32=True)
    read_idXML(args.input, args.output, where=args.where,
               no_run_group=args.no_run_group, filters=filters)


if __name__ == '__main__':
    main()
