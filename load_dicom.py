import pydicom
from dataclasses import dataclass
import pathlib
import numpy as np
import plotly.graph_objects as go


@dataclass
class dcm_lite:
    series: int
    imgindx: int
    pixelarr: np.ndarray
    header: dict


def create_header_dict(dcm_values):
    return


def create_dcm_array(datadirstr: str):
    dat = dict()
    datapath = pathlib.Path(datadirstr)

    for dcmf in datapath.iterdir():
        dcm = pydicom.dcmread(dcmf)
        header = create_header_dict(dcm.values())
        series_tag = ["0020", "0011"]
        imgindx_tag = ["0054", "1330"]
        key = dcm[series_tag].value
        if imgindx_tag in dcm:
            val = dcm[imgindx_tag].value
            dcm_obj = dcm_lite(key, val, dcm.pixel_array(), header)
            if key in dat:
                dat[key].append(dcm_obj)
            else:
                dat[key] = [dcm_obj]
        else:
            dat["nonseries"].append(
                dcm_lite(None, None, dcm.pixel_array(), header)
                )
    for series in dat:
        dat[series] = sorted(dat[series], key=lambda x: x.imgindx)

    return dat


def frame_args(duration):
    return {
        "frame": {"duration": duration},
        "mode": "immediate",
        "fromcurrent": True,
        "transition": {"duration": duration, "easing": "linear"},
    }


def make_fig(data):

    r, c = data[0].shape
    nb_frames = len(data)
    znum = data.shape[-1] * 0.1
    fig = go.Figure(
        frames=[
            go.Frame(
                data=go.Surface(
                    z=(znum - k * 0.1) * np.ones((r, c)),
                    surfacecolor=np.flipud(data[(len(data) - 1) - k]),
                    cmin=0,
                    cmax=25_000,
                ),
                name=str(k),
            )
            for k in range(nb_frames)
        ]
    )

    fig.add_trace(
        go.Surface(
            z=znum * np.ones((r, c)),
            surfacecolor=np.flipud(data[-1]),
            colorscale="gray_r",
            cmin=0,
            cmax=25_000,
            colorbar=dict(thickness=20, ticklen=4),
        )
    )

    sliders = [
        {
            "pad": {"b": 10, "t": 60},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [[f.name], frame_args(0)],
                    "label": str(k),
                    "method": "animate",
                }
                for k, f in enumerate(fig.frames)
            ],
        }
    ]

    # Layout
    fig.update_layout(
        title="Slices in volumetric data",
        width=600,
        height=600,
        scene=dict(
            zaxis=dict(autorange=True),
            aspectratio=dict(x=1, y=1, z=1),
        ),
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None, frame_args(50)],
                        "label": "&#9654;",  # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;",  # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
        ],
        sliders=sliders,
    )
    return fig
