import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from flask import Flask,session
import numpy as np
from SIR.advanced_SEIRS import simulateSEIRS,Progresses
from SIR import advanced_SEIRS
from scipy import interpolate
import secrets
import inspect


def crange(start, end, step=1.0):
    return np.arange(start, end + step / 2, step)


class Param:
    def __init__(self, *, name, min, max, primary_step=1.0, secondary_step=None, default=None, **kwargs):
        if False in (isnumeric := [type(i) in [int, float] for i in [min, max, primary_step, secondary_step]]):
            if secondary_step is not None:
                raise ValueError(f'{isnumeric.count(False)} of the parameters is not valid.')

        self.name = name
        self.min = min
        self.max = max
        self.primary_step = primary_step
        self.secondary_step = secondary_step if secondary_step is not None else primary_step
        self.default = default
        self.__dict__.update(kwargs)

    def array(self):
        return crange(self.min, self.max, self.secondary_step)

    def pmy_array(self):
        return crange(self.min, self.max, self.primary_step)
    def default_val(self):
        return self.default if self.default is not None else self.mid()
    def mid(self):
        return (x := self.array())[int(x.shape[0] / 2)]

server = Flask(__name__)  # real <'Flask'> object
app = dash.Dash(server=server)  # Dash object containing Flask object

app.title = '傳染病模擬(SERIS模型)'
server.config['SECRET_KEY']=secrets.token_bytes(16)


params = {
    'H_MAX': Param(
        min=0,
        max=500,
        primary_step=50,
        secondary_step=1,
        default=20,
        name='H_MAX'
    ),
    'beta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='beta'
    ),
    'sigma': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='sigma'
    ),
    'lambda_': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='lambda_'
    ),
    'a': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='a'
    ),
    'gamma_a': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='gamma_a'
    ),
    'h': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='h'
    ),
    'eta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='eta'
    ),
    'f_s': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='f_s'
    ),
    'mu_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='mu_s'
    ),
    'gamma_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='gamma_s'
    ),
    'f_h': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='f_h'
    ),
    'mu_h': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='mu_h'
    ),
    'gamma_h': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='gamma_h'
    ),
    'xi': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        default=0,
        name='x1'
    ),
    'Q_sigma': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_sigma'
    ),
    'Q_lambda_': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_lambda_'
    ),
    'Q_gamma_a': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_gamma_a'
    ),
    'Q_eta': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_eta'
    ),
    'Q_mu_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_mu_s'
    ),
    'Q_gamma_s': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='Q_gamma_s'
    ),
    'rho': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.1,
        name='rho'
    ),
    'theta_E': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.01,
        name='theta_E'
    ),
    'psi_E': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='psi_E'
    ),
    'theta_I_pre': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.01,
        name= 'theta_I_pre'
    ),
    'psi_I_pre': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='psi_I_pre'
    ),
    'theta_I_asym': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.01,
        name='theta_I_asym'
    ),
    'psi_I_asym': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='psi_I_asym'
    ),
    'theta_I_sym': Param(
        min=0,
        max=3,
        primary_step=1,
        secondary_step=0.01,
        name='theta_I_sym'
    ),
    'psi_I_sym': Param(
        min=0,
        max=1,
        primary_step=0.1,
        name='psi_I_sym'
    ),
    'day': Param(
        min=0,
        max=365,
        primary_step=20,
        secondary_step=1,
        name='day'
    )
}
S = 150
I_sym = 135
F = 15



def create_layout():
    layout = html.Div([
        html.H2(
            id='title',
            children='SEIRS Simulation'
        ),
        html.H3(
            id='warning-label',
            children="頁面下方有說明"
        ),
        html.Div([
            dcc.Dropdown(
                id='progress-dropdown',
                options=[
                    {'label': f'{progress[0]}: '+progress[1].__doc__.replace('\n','  '), 'value': (val:=progress[0])}
                    for progress in inspect.getmembers(Progresses)\
                    if not (progress[0].startswith('__') and progress[0].endswith('__')) and progress[0].find('progress_')!=-1
                ],
                value=val,
                placeholder="Select a Progress",
            ),
            html.Div(id='dd-output-container',style={'margin-top':'10px'})
        ],style={'margin':'10px 0px'}),
        dcc.Graph(
            id='seirs',
            animate=True,
            style = {'margin': '20px 0px'}
        ),
        *[
            html.Div([
                dcc.Slider(
                    id=key,
                    min=param.min,
                    max=param.max,
                    value=param.default_val(),
                    step=param.secondary_step,
                    marks={(round(val,10) if not val==int(val) else int(val)): {'label': str(round(val,10))} for val in param.pmy_array()},
                    persistence=True
                ),
                html.Div(id=f'{key}-label',style={'margin-bottom':30,'text-align':"center"})
            ]) for key,param in params.items()
        ],
        html.Div([
            html.H4(children='說明'),
            dcc.Markdown(advanced_SEIRS.__doc__,
                         style={'white-space':'pre','overflow-x':'scroll'})
        ],style={'padding':"0px 10vw"})

    ])
    return layout

app.layout = create_layout()

@app.callback(
    dash.dependencies.Output('seirs', 'figure'),
    [
        *[dash.dependencies.Input(key, 'value') for key in params.keys()],
        dash.dependencies.Input('progress-dropdown','value')
    ]
)
def update_figure(*vals):
    traces = []
    print(vals)
    progress=vals[-1]
    print(f'using progress "{progress}" to similate!')
    vals=vals[:-1]
    keys_n_vals = {key: val for key, val in zip(params.keys(), vals)}
    session.update(keys_n_vals)
    SEIRS, params_used = simulateSEIRS(**keys_n_vals,S=S,I_sym=I_sym,F=F,progress_func=dict(inspect.getmembers(Progresses))[progress])
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
        'E':SEIRS[:, params_used['compartments'].index('E')],
        'I':np.sum(SEIRS[:,[i for i in range(types_amount) if True in list(map(lambda symbol:symbol in params_used['compartments'][i], ['I'] )) ]], axis=1),
        'H':np.sum(SEIRS[:,[i for i in range(types_amount) if True in list(map(lambda symbol:symbol in params_used['compartments'][i], ['H'] )) ]], axis=1),
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
            xaxis={'title': 'days', 'range': [0, params_used['day']]},
            yaxis={'title': 'population(people)', 'range': [0, np.max(SEIRS[:, params_used['compartments'].index('total')]) * 1.2]},
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

@app.callback(
    dash.dependencies.Output('dd-output-container', 'children'),
    [dash.dependencies.Input('progress-dropdown', 'value')])
def update_output(value):
    return f'You have selected "{value}"'

if __name__ == '__main__':
    app.run_server(debug=True)
