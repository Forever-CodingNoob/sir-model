from pylab import *
from scipy import interpolate
import numpy as np
import time
import math

# Parameters
Beta = 0.6 #beta
Gamma = 0.1 #gamma
Mu = 0.1 #疾病致死率
Alpha = 1/7 #1/潛伏期
N = 300
day = 30
interval = 0.2
I0 = 135
E0 = 0
R0 = 15
D0 = 0
# birth = 0 # 出生率 = 新出生人口數 / 原總人數 = (新總人數-原總人數) / 原總人數 (若不考慮死亡率)
# death = 0 # 死亡率


names = ['S','E','I','R','D']

def runSEIRD(S, E, I, R, D, dt, t):
    global N
    # dBorndt = (1/365)*(N-D)*((birth+1)**0)*math.log(birth+1)
    dSdt = - Beta * S * I / N
    dRdt = Gamma * I
    dDdt = Mu * I
    dEdt = - Alpha*E
    dIdt = - dEdt - dRdt - dDdt
    S += (dSdt) * dt
    E += (-dSdt+dEdt) * dt
    I += dIdt * dt
    R += dRdt * dt
    D += dDdt * dt

    # S,I,R,D=map(int,[S,I,R,D]) # turn S,I,R,D into integer

    def fix(n):
        if n<0:
            return 0
        return n
    S,E,I,R,D=map(fix,[S,E,I,R,D])

    if I<=0:
        print('disease gone!')

    # N =S+E+I+R+D
    return S, E, I, R, D


def integrateSEIRD(SEIRD, total_time, interval):
    t = 0
    for i in range(int(total_time / interval)):
        new_SEIRD = np.array(runSEIRD(*SEIRD[-1], interval, t))
        # print(new_SEIRD)
        # print(SEIRD)
        SEIRD = np.vstack((SEIRD, new_SEIRD))

        t += interval

    return SEIRD


if __name__ == '__main__':
    now = time.time()

    '''運行SIR模擬'''
    E = E0
    I = I0
    R = R0
    D = D0
    S = N - I0 - R0 - E0 - D0


    SEIRD = np.array([[S, E, I, R, D]], dtype=float)

    SEIRD = integrateSEIRD(SEIRD, day, interval)
    print(SEIRD[-1], N)

    '''畫圖'''
    length = np.size(SEIRD, axis=0)
    types_amount = np.size(SEIRD, axis=1)
    t = arange(0, length, 1) * (interval)  # [0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻
    # print(id(LS))
    # print(t)
    # print(LR.size)
    # print(t.size)

    # 運用內差法
    # 繪製曲線
    dt = 0.05
    tnew = arange(0, int(day / dt) + 1, 1) * dt  # [0,0+dt,0+dt*2,...,day]

    lines = [plot(tnew, interpolate.InterpolatedUnivariateSpline(t, SEIRD[:,i])(tnew), label=names[i])[0] for i in range(types_amount)]
    print(lines)

    # 圖例
    legend(handles=lines,
           shadow=True, loc=(0.85, 0.4))  # handle

    # 標示
    text(16.5, 290, 'Beta=%g' % (Beta))
    text(16.5, 270, 'Gamma=%g' % (Gamma))
    text(16.5, 250, 'Mu=%g' % (Mu))
    text(16.5, 230, 'Alpha=%g' % (Alpha))
    text(16.5, 210, 'Infected:%d' % (I0))

    v = [0, day, 0, N]  # [x刻度start,x刻度end,y刻度start,y刻度end]
    axis(v)
    xlabel('time (days)')
    ylabel('Population(people)')
    title('SEIRD Model')
    grid(True)
    print('execution time: %fs' % (time.time() - now))
    show()