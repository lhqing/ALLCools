"""
- Input: a tsv file that contain layout settings and file information
- Validate all params
- Make a session object
- Output: IGV session XML
"""

from itertools import product
from .Session import Session
import pandas as pd
from .defaults import *
import xml.etree.ElementTree as ET
import pathlib
from urllib.request import urlopen


def hex_to_rgb_text(color_code):
    rgb = tuple(int(color_code.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4))
    return ','.join(map(str, rgb))


def sort_df_by_order_dict(df, order_dict):
    df = df.set_index(list(order_dict.keys())).reindex(
        list(product(*order_dict.values()))).reset_index(
    )  # sort by product of all order list
    return df


def make_session(track_table, output_xml_path,
                 hue=None, palette=None,
                 range_key_col=None, data_range_dict=None,
                 genome='mm10', locus='Gad1', open_in_igv=True):
    track_df = pd.read_csv(track_table)

    not_in = []
    for col in TABLE_REQUIRED_FIELDS:
        if col not in track_df.columns:
            not_in.append(col)
    if len(not_in) != 0:
        raise KeyError(f'Missing {not_in} col in track_df')

    if hue is not None:
        if palette is not None:
            track_df['color'] = track_df[hue].apply(lambda i: hex_to_rgb_text(palette[i]))
        else:
            raise ValueError('Provide palette together with hue.')

    # init igv session
    session = Session(
        genome=genome,
        hasGeneTrack='true',
        hasSequenceTrack='true',
        locus=locus,
        version='8',
        refseq_track=True)

    # add track
    for _, row in track_df.iterrows():
        panel = row.pop('panel')
        file_format = row.pop('file_format')
        if file_format.lower() in FEATURE_TRACK_FORMATS:
            track_type = 'Feature'
            default_dict = FEATURE_TRACK_DEFAULT
        else:
            track_type = 'Data'
            default_dict = DATA_TRACK_DEFAULT

        track_kws = {'name': row['modality'] + ' - ' + row['name'],
                     'id': row['path']}
        for k, v in default_dict.items():
            if k in row:
                track_kws[k] = str(row[k])
            else:
                track_kws[k] = v

        # data range only apply to data track
        data_range_kws = {}
        for k, v in DATARANGE_DEFAULT.items():
            if k in row:
                data_range_kws[k] = str(row[k])
            else:
                data_range_kws[k] = v
            if (range_key_col is not None) and (row[range_key_col] in data_range_dict):
                data_range_kws['minimum'], data_range_kws['maximum'] = map(str, data_range_dict[row[range_key_col]])

        if track_type == 'Feature':
            session.add_feature_track(track_kws=track_kws, panel=panel)
        else:
            session.add_data_track(track_kws=track_kws, data_range_kws=data_range_kws, panel=panel)

    with open(output_xml_path, 'w') as f:
        f.write(ET.tostring(session, encoding='utf8').decode('utf8'))

    if open_in_igv:
        output_xml_path = pathlib.Path(output_xml_path).absolute()
        response = urlopen(f"http://localhost:60151/load?file={output_xml_path}")
        response.read(100)
    return
