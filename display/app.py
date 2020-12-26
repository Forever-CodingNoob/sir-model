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


params = {
    'beta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    )
}



app.layout = html.Div([
    html.H2(children='SEIRS'),
    dcc.Graph(id='seirs',
              animate=True
              ),
    dcc.Slider(
        id='beta',
        min=params['beta'].min,
        max=params['beta'].max,
        value=params['beta'].mid(),
        step=params['beta'].secondary_step,
        marks={(round(val,10) if not val==int(val) else int(val)): {'label': str(round(val,10))} for val in params['beta'].array()}
    )
])


@app.callback(
    dash.dependencies.Output('seirs', 'figure'),
    [dash.dependencies.Input('beta', 'value')])
def update_figure(beta):
    traces = []
    print('beta:', beta)
    SEIRS, params = simulateSEIRS(beta=beta,day=30)
    print(SEIRS, params, sep='\n')

    length = np.size(SEIRS, axis=0)
    types_amount = np.size(SEIRS, axis=1)
    t = np.arange(0, length, 1) * (params['interval'])  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻

    # 運用內差法
    # 繪製曲線
    dt = 0.05
    tnew = np.arange(0, int(params['day'] / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]
    # print(tnew)

    traces = [go.Scatter(
        x=tnew,
        y=interpolate.InterpolatedUnivariateSpline(t, SEIRS[:, i])(tnew),
        text=params['compartments'][i],
        mode='lines',
        opacity=0.7,
        name=params['compartments'][i]
    ) for i in range(types_amount)]

    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'title': 'day', 'range': [0,params['day']]},
            yaxis={'title': 'population', 'range': [0,np.max(SEIRS[:, -1]) * 1.2]},
            margin={'l': 100, 'b': 40, 't': 10, 'r': 10},
            legend=go.layout.Legend(),
            hovermode='closest'
        )
    }


if __name__ == '__main__':
    app.run_server(debug=True)
