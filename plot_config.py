from pathlib import Path
from toolz import pipe

import altair as alt
import matplotlib as mpl

JSON_DATA_PATH = Path('altair-data/')


def json_dir(data, data_dir: Path = JSON_DATA_PATH):
    data_dir.mkdir(exist_ok=True)
    return pipe(
        data, alt.to_json(filename=str(data_dir / '{prefix}-{hash}.{extension}'))
    )


def init_json_transformer():
    alt.data_transformers.disable_max_rows()
    alt.data_transformers.register('json_dir', json_dir)
    alt.data_transformers.enable('json_dir', data_dir=JSON_DATA_PATH)


def altair_figure_theme(*args, **kwargs):
    return {'font': 'Arial', 'fontSize': 10}


def config_svg():
    new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
    mpl.rcParams.update(new_rc_params)
