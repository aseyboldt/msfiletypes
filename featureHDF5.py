import xml.etree.ElementTree as ET
import tables as tb
import argparse
import contextlib
import numpy as np


class Features(tb.IsDescription):
    """ Corresponds to a feature in featureXML. """
    id = tb.StringCol(40, pos=0)
    x_pos = tb.Float64Col(pos=1)
    y_pos = tb.Float64Col(pos=2)
    intensity = tb.Float64Col(pos=3)
    overallquality = tb.Float64Col(pos=4)
    charge = tb.Int32Col(pos=5)
    score_fit = tb.Float64Col(pos=6)
    score_correlation = tb.Float64Col(pos=7)
    fwhm = tb.Float64Col(pos=8)
    spectrum_index = tb.Int64Col(pos=9)


def write_features(iterparser, start_tag, table):

    def read_feature(elem, buffer):
        buffer['id'] = elem.get('id')
        buffer['x_pos'] = float(elem.find("position[@dim='0']").text)
        buffer['y_pos'] = float(elem.find("position[@dim='1']").text)
        buffer['intensity'] = float(elem.find("intensity").text)
        buffer['overallquality'] = float(elem.find("overallquality").text)
        buffer['charge'] = int(elem.find("charge").text)
        buffer['score_fit'] = float(
            elem.find("userParam[@name='score_fit']").get('value')
        )
        buffer['score_correlation'] = float(
            elem.find("userParam[@name='score_correlation']").get('value')
        )
        buffer['fwhm'] = float(
            elem.find("userParam[@name='FWHM']").get('value')
        )
        buffer['spectrum_index'] = int(
            elem.find("userParam[@name='spectrum_index']").get('value')
        )

    buffer_size = 1024
    buffer_array = np.empty((buffer_size, ), table.dtype)

    buffer_index = 0
    for event, elem in iterparser:
        if event == 'start':
            continue
        if elem.tag == 'featureList':
            table.append(buffer_array[:buffer_index])
            start_tag.clear()
            return
        elif elem.tag == 'feature':
            read_feature(elem, buffer_array[buffer_index])
            buffer_index += 1
            if buffer_index == buffer_size:
                table.append(buffer_array)
                buffer_index = 0
            start_tag.clear()


def write_feature_list(iterparser, start_tag, hdf_file,
                       where="/", filters=None):
    """ Store the data from one featureList in a hdf5 group.

    Parameters
    ----------
    search_params: iterator like ElementTree.iterparse
        This must be in a state after ('start', 'FeatureList') has been
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
    if not start_tag.tag == 'featureList':
        raise ValueError('Not a featureList tag: %s' % start_tag.tag)

    if where not in hdf_file:
        where = hdf_file.create_group('/', where.lstrip('/'))

    feature_count = None
    if 'count' in start_tag.attrib:
        feature_count = int(start_tag.attrib['count'])
    features = hdf_file.create_table(
        where, 'features', description=Features, filters=filters,
        expectedrows=feature_count
    )

    write_features(iterparser, start_tag, features)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an OpenMS featureXML file to a hdf5 format"
    )

    parser.add_argument('input', help='input idXML file', type=str)
    parser.add_argument('output', help='Output path')
    parser.add_argument('--compress', choices=['blosc', 'zlib', 'lzo'],
                        default='zlib')
    parser.add_argument('--compression-level', choices=list(range(1, 10)),
                        default=1, type=int)
    parser.add_argument('--where', help='hdf5 destination inside output',
                        default='/')
    return parser.parse_args()


def main():
    args = parse_args()
    infile = ET.iterparse(args.input, ['start', 'end'])

    filters = tb.Filters(complevel=args.compression_level,
                         complib=args.compress, fletcher32=True)
    f = tb.open_file(args.output, 'a')
    with contextlib.closing(f):
        for event, elem in infile:
            if event == 'start' and elem.tag == 'featureList':
                write_feature_list(
                    infile, elem, f, args.where, filters=filters
                )

if __name__ == '__main__':
    main()
