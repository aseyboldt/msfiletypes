import xml.etree.ElementTree as ET
import tables as tb
import argparse
import contextlib
import numpy as np
import h5py
import base64


grouping_accessions = [
    'MS:1000523',  # float64
    'MS:1000521',  # float32
    'MS:1000128',  # profile spectrum
    'MS:1000511',  # ms level, has value
    'MS:1000579',  # MS1 spectrum
]


def write_data_array(elem, where, name):
    data = elem.find('{http://psi.hupo.org/ms/mzml}binary')
    type = elem.find('{http://psi.hupo.org/ms/mzml}cvParam[@accession="MS:1000523"]')
    if type is not None:
        dtype = 'd'
    elif elem.find('{http://psi.hupo.org/ms/mzml}cvParam[@accession="MS:1000521"]') is not None:
        dtype = 'f'
    else:
        raise ValueError("unknown data type: %s" % type.attr['name'])
    decoded = base64.b64decode(data.text)
    array = np.frombuffer(decoded, dtype=dtype)
    where.create_dataset(name, compression="gzip", shuffle=True, data=array)


def convert_mzml(infile, storage, where='/'):
    """ Read mzml and store the arrays in hdf storage.

    Parameters
    ----------
    infile : path or file
        Input mzml file to read
    storage : path or h5py.File instance
        The HDF5 file to store the tables in. If a path is specified the
        file will be opened in append mode.
    where : str or pytables.Group
        The path in the HDF5 file to use
    """
    def read_groups(storage_file):
        if where not in storage_file:
            storage_file.create_group(where)
        array_id = 0
        for event, elem in infile:
            if event == 'end' and elem.tag == '{http://psi.hupo.org/ms/mzml}binaryDataArray':
                write_data_array(
                    elem,
                    where=storage_file[where],
                    name='data_%08i' % array_id
                )
                array_id += 1
            #root.clear()

    infile = ET.iterparse(infile, ['start', 'end'])
    _, root = next(infile)
    if isinstance(storage, h5py.File):
        read_groups(storage)
    else:
        f = h5py.File(storage, 'a')
        with contextlib.closing(f):
            read_groups(f)
Quoting Hannes Röst <hannesroest@gmail.com>:

> Hi Adrian
>
> Wie versprochen noch der Link zu der Doku fuer den file access
> http://ftp.mi.fu-berlin.de/OpenMS/documentation/html/tutorial_format.html
> und als Anhang noch das chapter meiner thesis das ich gerade dazu
> schreibe. Bin gespannt von dir zu hoeren wenn es Neuigkeiten zu hdf5
> gibt.
>
> Beste Gruesse
>
> Hannes


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an mzML file to a hdf5 format"
    )

    parser.add_argument('input', help='input mzML file', type=str)
    parser.add_argument('output', help='Output path')
    parser.add_argument('--compress', choices=['blosc', 'zlib', 'lzo', 'lzf'],
                        default='zlib')
    parser.add_argument('--compression-level', choices=list(range(1, 10)),
                        default=1, type=int)
    parser.add_argument('--where', help='hdf5 destination inside output',
                        default='/')
    return parser.parse_args()


def main():
    args = parse_args()
    convert_mzml(args.input, args.output, where=args.where)


if __name__ == '__main__':
    main()
