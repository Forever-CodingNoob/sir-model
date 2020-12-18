from pylab import *
from scipy import interpolate
import numpy as np
import time

# Parameters
Beta = 0.1
Gamma = 0.05
N = 1000
day = 365
interval = 0.2
I0 = 5


def SIR(S, I, R, dt):

    dS = - Beta*S*I/N
    dR = Gamma*I
    dI = - dS - dR
    S += dS*dt
    I += dI*dt
    R += dR*dt
    if S < 0:
        S = 0
    if I > N:
        I = N
    if R > N:
        R = N

    return S, I, R
def integrateSIR(LS,LI,LR,total_time,interval):
    t=0
    for i in range(int(total_time/interval)):
        s,i,r=SIR(LS[-1],LI[-1],LR[-1],interval)

        # print(id(LS))
        LS = np.append(LS,s) #`append` creates a new ndarray different from the former `LS` ndarray then assign the `LS` variable name to the newly-built array
        LI = np.append(LI,i) # hence the `LS` here is different (having a distinct id) from the one in main function
        LR = np.append(LR,r) # it's necessary to return the 3 arrays in this function, otherwise the ones in main function are not even modified or reassigned
        # print(t)

        t+=interval

    return LS,LI,LR



if __name__ == '__main__':
    now=time.time()

    '''運行SIR模擬'''
    I = I0
    S = N - I0
    R = 0

    LS = np.array([S],dtype=float)
    LI = np.array([I],dtype=float)
    LR = np.array([R],dtype=float)

    LS,LI,LR=integrateSIR(LS, LI, LR, day, interval)

    '''畫圖'''
    t = arange(0, LS.size, 1)*(interval) #[0,0+interval,0+interval*2,...,day]，即每筆資料對應之時刻
    # print(id(LS))
    # print(t)
    # print(LR.size)
    # print(t.size)

    # 運用內差法
    ls = interpolate.InterpolatedUnivariateSpline(t, LS)
    li = interpolate.InterpolatedUnivariateSpline(t, LI)
    lr = interpolate.InterpolatedUnivariateSpline(t, LR)

    # 繪製曲線
    dt=0.05
    tnew = arange(0, int(day/dt)+1,1)*dt #[0,0+dt,0+dt*2,...,day]
    lsnew = ls(tnew)
    linew = li(tnew)
    lrnew = lr(tnew)
    line1, = plot(tnew, lsnew, label='S')
    line2, = plot(tnew, linew, label='I')
    line3, = plot(tnew, lrnew, label='R')

    # 圖例
    legend(handles=[line1, line2, line3],
           shadow=True, loc=(0.85, 0.4))  # handle

    # 描上資料點
    # line11, = plot(t, LS, '.')
    # line22, = plot(t, LI, '.')
    # line33, = plot(t, LR, '.')

    # 標示
    text(16.5, 240, 'Beta=%g' % (Beta))
    text(16.5, 190, 'Gamma=%g' % (Gamma))
    text(16.5, 150, 'Infected:%d' % (I0))

    v = [0, day, 0, 1000] #[x刻度start,x刻度end,y刻度start,y刻度end]
    axis(v)
    xlabel('time (days)')
    ylabel('Population(people)')
    title('SIR Model')
    grid(True)
    print('execution time: %fs' %(time.time()-now))
    show()