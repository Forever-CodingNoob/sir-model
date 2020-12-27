import dash
import dash_core_components as dcc
import dash_html_components as html
import math
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import typing
from SIR.advanced_SEIRS import simulateSEIRS
from scipy import interpolate


def crange(start, end, step=1.0):
    return np.arange(start, end + step / 2, step)


class Param:
    def __init__(self, *, min, max, primary_step=1.0, secondary_step=None, **kwargs):
        if False in (isnumeric := [type(i) in [int, float] for i in [min, max, primary_step, secondary_step]]):
            raise ValueError(f'{isnumeric.count(False)} of the parameters is not valid.')

        self.min = min
        self.max = max
        self.primary_step = primary_step
        self.secondary_step = secondary_step if secondary_step is not None else primary_step
        self.__dict__.update(kwargs)

    def array(self):
        return crange(self.min, self.max, self.secondary_step)

    def pmy_array(self):
        return crange(self.min, self.max, self.primary_step)

    def mid(self):
        return (x := self.array())[int(x.shape[0] / 2)]


app = dash.Dash()
server = app.server


params = {
    'beta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'I_sym': Param(
        min=0,
        max=150,
        primary_step=10,
        secondary_step=1
    ),
    'R': Param(
        min=0,
        max=150,
        primary_step=10,
        secondary_step=1
    )
}



app.layout = html.Div([
    html.H2(children='SEIRS'),
    dcc.Graph(id='seirs',
              animate=True
              ),
    *[
        html.Div([
            dcc.Slider(
                id=key,
                min=param.min,
                max=param.max,
                value=param.mid(),
                step=param.secondary_step,
                marks={(round(val,10) if not val==int(val) else int(val)): {'label': str(round(val,10))} for val in param.pmy_array()}
            ),
            html.Div(id=f'{key}-label',style={'margin-bottom':20})
        ]) for key,param in params.items()
    ]

])


@app.callback(
    dash.dependencies.Output('seirs', 'figure'),
    [dash.dependencies.Input(key, 'value') for key in params.keys()])
def update_figure(*vals):
    traces = []
    print(vals)
    SEIRS, params_used = simulateSEIRS(**{key: val for key, val in zip(params.keys(), vals)}, day=30)
    print(SEIRS, params_used, sep='\n')

    length = np.size(SEIRS, axis=0)
    types_amount = np.size(SEIRS, axis=1)
    t = np.arange(0, length, 1) * (params_used['interval'])  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻

    # 運用內差法
    # 繪製曲線
    dt = 0.05
    tnew = np.arange(0, int(params_used['day'] / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]
    # print(tnew)

    traces = [go.Scatter(
        x=tnew,
        y=interpolate.InterpolatedUnivariateSpline(t, SEIRS[:, i])(tnew),
        mode='lines',
        opacity=0.7,
        name=params_used['compartments'][i]
    ) for i in range(types_amount)]

    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'title': 'day', 'range': [0, params_used['day']]},
            yaxis={'title': 'population', 'range': [0, np.max(SEIRS[:, -1]) * 1.2]},
            margin={'l': 100, 'b': 40, 't': 10, 'r': 10},
            legend=go.layout.Legend(),
            hovermode='closest'
        )
    }

def gen_updLabel_func(keyname):
    def upd_label(val):
        return keyname+" = "+str(val)
    return upd_label
for key,param in params.items():
    # upd_label.__closure__[0].cell_contents = key

    app.callback(
        dash.dependencies.Output(f'{key}-label', 'children'),
        [dash.dependencies.Input(key, 'value')]
    )(gen_updLabel_func(keyname=key))



if __name__ == '__main__':
    app.run_server(debug=True)
