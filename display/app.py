import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import numpy as np
from SIR.advanced_SEIRS import simulateSEIRS,progress_01,progress
from scipy import interpolate


def crange(start, end, step=1.0):
    return np.arange(start, end + step / 2, step)


class Param:
    def __init__(self, *, min, max, primary_step=1.0, secondary_step=None, **kwargs):
        if False in (isnumeric := [type(i) in [int, float] for i in [min, max, primary_step, secondary_step]]):
            if secondary_step is not None:
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
    'sigma': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'lambda_': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'a': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'gamma_a': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'h': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'eta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'f_s': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'mu_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'gamma_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'f_h': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'mu_h': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'gamma_h': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'xi': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_sigma': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_lambda_': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_gamma_a': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_eta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_mu_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'Q_gamma_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'rho': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'theta_E': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'psi_E': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'theta_I_pre': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'psi_I_pre': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'theta_I_asym': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'psi_I_asym': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'theta_I_sym': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1
    ),
    'psi_I_sym': Param(
        min=0,
        max=1,
        primary_step=0.1
    ),
    'day': Param(
        min=0,
        max=365,
        primary_step=20,
        secondary_step=1
    )
}
S = 150
I_sym = 135
F = 15



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
    SEIRS, params_used = simulateSEIRS(**{key: val for key, val in zip(params.keys(), vals)},S=S,I_sym=I_sym,F=F,progress_func=progress_01)
    print(SEIRS, params_used, sep='\n')

    length = np.size(SEIRS, axis=0)
    types_amount = np.size(SEIRS, axis=1)
    t = np.arange(0, length, 1) * (params_used['interval'])  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻

    '''
    運用內差法繪製曲線
    dt = 0.05
    tnew = np.arange(0, int(params_used['day'] / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]
    print(tnew)


    traces = [go.Scatter(
        x=tnew,
        y=interpolate.InterpolatedUnivariateSpline(t, SEIRS[:, i])(tnew),
        mode='lines',
        opacity=0.7,
        name=params_used['compartments'][i]
    ) for i in range(types_amount)]
    '''
    combined_SEIRS={
        'S':SEIRS[:, params_used['compartments'].index('S')],
        'I':np.sum(SEIRS[:,[i for i in range(types_amount) if True in list(map(lambda symbol:symbol in params_used['compartments'][i], ['E','I','H'] )) ]], axis=1),
        'R':SEIRS[:, params_used['compartments'].index('R')],
        'F':SEIRS[:, params_used['compartments'].index('F')],
        'Q':np.sum(SEIRS[:,[i for i in range(types_amount) if 'Q' in params_used['compartments'][i]]], axis=1),
        'total(N)':SEIRS[:, params_used['compartments'].index('total')]
    }
    # traces = [go.Scatter(
    #     x=t,
    #     y=SEIRS[:, i],
    #     mode='lines',
    #     opacity=0.7,
    #     name=params_used['compartments'][i]
    # ) for i in range(types_amount)]
    traces = [go.Scatter(
        x=t,
        y=data,
        mode='lines',
        opacity=0.7,
        name=name
    ) for name,data in combined_SEIRS.items()]

    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'title': 'day', 'range': [0, params_used['day']]},
            yaxis={'title': 'population', 'range': [0, np.max(SEIRS[:, params_used['compartments'].index('total')]) * 1.2]},
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
