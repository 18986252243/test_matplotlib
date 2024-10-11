# -*- coding: utf-8 -*-

from math import sqrt, fabs
from typing import Dict, List, Tuple
import matplotlib as mpl
import matplotlib.pyplot as plt

meas_data_list:List = [
    {'pos_id':'TRP0',   'd_A':5.06,     'd_B':13.60,    'd_C':18.87,    'H':1.62,       'Ref_A':'pos3_3',     'Ref_B':'pos10_8',     'Ref_C':'pos3_18',     'x_cor':None,   'y_cor':None,   'err':None,   'Sel':None},
    {'pos_id':'TRP1',   'd_A':16.68,    'd_B':15.82,    'd_C':4.350,    'H':1.60,       'Ref_A':'pos3_3',     'Ref_B':'pos10_8',     'Ref_C':'pos3_18',     'x_cor':None,   'y_cor':None,   'err':None,   'Sel':None},
    {'pos_id':'TRP2',   'd_A':19.09,    'd_B':14.48,    'd_C':26.35,    'H':1.57,       'Ref_A':'pos3_3',     'Ref_B':'pos10_8',     'Ref_C':'pos3_18',     'x_cor':None,   'y_cor':None,   'err':None,   'Sel':None},
    {'pos_id':'TRP3',   'd_A':24.95,    'd_B':16.47,    'd_C':18.42,    'H':1.57,       'Ref_A':'pos3_3',     'Ref_B':'pos10_8',     'Ref_C':'pos3_18',     'x_cor':None,   'y_cor':None,   'err':None,   'Sel':None},
]

ref_point_list:List = [
    {'pos_id':'pos3_3',   'x_cor':4.67,     'y_cor':4.63},
    {'pos_id':'pos10_8',  'x_cor':11.72,    'y_cor':9.52},
    {'pos_id':'pos3_18',  'x_cor':4.83,     'y_cor':19.62},
]

test_field_size:Dict = {'x_size':35.2, 'y_size':22}

def get_C(meas_data_dict:Dict, ref_A_dict:Dict, ref_B_dict:Dict) -> float:
    d_A = meas_data_dict['d_A']
    d_B = meas_data_dict['d_B']

    x_A = ref_A_dict['x_cor']
    y_A = ref_A_dict['y_cor']
    x_B = ref_B_dict['x_cor']
    y_B = ref_B_dict['y_cor']
     
    C = ((d_A**2-d_B**2) - (x_A**2-x_B**2 + y_A**2-y_B**2))/(2*(x_B-x_A))
    
    return C

def get_K0(meas_data_dict:Dict, ref_A_dict:Dict, C:float) -> float:
    d_A = meas_data_dict['d_A']
    x_A = ref_A_dict['x_cor']
    y_A = ref_A_dict['y_cor']
    
    K0 = C**2 - 2*x_A*C + x_A**2 + y_A**2 - d_A**2
    
    return K0

def get_K1(meas_data_dict:Dict, ref_A_dict:Dict, ref_B_dict:Dict, C:float) -> float:
    d_A = meas_data_dict['d_A']
    d_B = meas_data_dict['d_B']

    x_A = ref_A_dict['x_cor']
    y_A = ref_A_dict['y_cor']
    x_B = ref_B_dict['x_cor']
    y_B = ref_B_dict['y_cor']
    
    K1 = 2*x_A*(y_B-y_A)/(x_B-x_A) - 2*C*(y_B-y_A)/(x_B-x_A) - 2*y_A
    
    return K1

def get_K2(ref_A_dict:Dict, ref_B_dict:Dict) -> float:
    x_A = ref_A_dict['x_cor']
    y_A = ref_A_dict['y_cor']
    x_B = ref_B_dict['x_cor']
    y_B = ref_B_dict['y_cor']
    
    K2 = 1 + ((y_B-y_A)**2)/((x_B-x_A)**2)
    
    return K2

def get_Y(meas_data_dict:Dict, ref_A_dict:Dict, ref_B_dict:Dict) -> List:
    Y = [None]*2
    
    C = get_C(meas_data_dict, ref_A_dict, ref_B_dict)
    K0 = get_K0(meas_data_dict, ref_A_dict, C)
    K1 = get_K1(meas_data_dict, ref_A_dict, ref_B_dict, C)
    K2 = get_K2(ref_A_dict, ref_B_dict)
    
    Y[0] = (0-K1-sqrt(K1**2-4*K2*K0))/(2*K2)
    Y[1] = (0-K1+sqrt(K1**2-4*K2*K0))/(2*K2)
    
    return Y

def get_X(meas_data_dict:Dict, ref_A_dict:Dict, ref_B_dict:Dict) -> List:
    X = [None]*2
    
    x_A = ref_A_dict['x_cor']
    y_A = ref_A_dict['y_cor']
    x_B = ref_B_dict['x_cor']
    y_B = ref_B_dict['y_cor']    
    C = get_C(meas_data_dict, ref_A_dict, ref_B_dict)
    Y = get_Y(meas_data_dict, ref_A_dict, ref_B_dict)
    
    X[0] = C - ((y_B-y_A)/(x_B-x_A))*Y[0]
    X[1] = C - ((y_B-y_A)/(x_B-x_A))*Y[1]
    
    return X

def get_Sel(meas_data_dict:Dict, ref_C_dict:Dict) -> Tuple:
    d0 = sqrt((meas_data_dict['x_cor'][0]-ref_C_dict['x_cor'])**2 + (meas_data_dict['y_cor'][0]-ref_C_dict['y_cor'])**2)
    d1 = sqrt((meas_data_dict['x_cor'][1]-ref_C_dict['x_cor'])**2 + (meas_data_dict['y_cor'][1]-ref_C_dict['y_cor'])**2)
    
    err0 = fabs(d0 - meas_p['d_C'])
    err1 = fabs(d1 - meas_p['d_C'])
    
    if err0 <= err1:
        Sel = 0
        err = d0 - meas_p['d_C']
    else:
        Sel = 1
        err = d1 - meas_p['d_C']
        
    return (Sel, err)

def print_table() -> None:
    global meas_data_list
    
    print('\tX[0]\tY[0]\tX[1]\tY[1]\tZ\tErr')
    for meas_p in meas_data_list:
        print('{}\t'.format(meas_p['pos_id']), end="")
        if meas_p['Sel'] == 0:
            print('*{:.3f}\t*{:.3f}\t{:.3f}\t{:.3f}\t{}\t{}'.format(meas_p['x_cor'][0], meas_p['y_cor'][0], meas_p['x_cor'][1], meas_p['y_cor'][1], meas_p['H'], meas_p['err']))
        else:
            print('{:.3f}\t{:.3f}\t*{:.3f}\t*{:.3f}\t{}\t{}'.format(meas_p['x_cor'][0], meas_p['y_cor'][0], meas_p['x_cor'][1], meas_p['y_cor'][1], meas_p['H'], meas_p['err']))
                
    return

def plot_ref_p(graph_handler:plt.Axes, ref_points:List) -> None:
    for ref_p in ref_points:
        graph_handler.scatter(ref_p['x_cor'], ref_p['y_cor'], label=ref_p['pos_id'], c='b', marker='s')
        plt.text(ref_p['x_cor'], ref_p['y_cor']+0.4, '{}\n({:.2f}, {:.2f})'.format(ref_p['pos_id'], ref_p['x_cor'], ref_p['y_cor']), ha='center')
    return    

def plot_meas_p(graph_handler:plt.Axes, meas_data_list:List) -> None:
    for meas_p in meas_data_list:
        graph_handler.scatter(meas_p['x_cor'][meas_p['Sel']], meas_p['y_cor'][meas_p['Sel']], label=meas_p['pos_id'])
        plt.text(meas_p['x_cor'][meas_p['Sel']], meas_p['y_cor'][meas_p['Sel']]+0.4, '{}\n({:.2f}, {:.2f})'.format(meas_p['pos_id'], meas_p['x_cor'][meas_p['Sel']], meas_p['y_cor'][meas_p['Sel']]), ha='center')
    return

# 计算坐标
for meas_p in meas_data_list:    
    for ref_p in ref_point_list:
        if ref_p['pos_id'] == meas_p['Ref_A']:
            ref_A = ref_p
        elif ref_p['pos_id'] == meas_p['Ref_B']:
            ref_B = ref_p
        elif ref_p['pos_id'] == meas_p['Ref_C']:
            ref_C = ref_p            
            
    meas_p['y_cor'] = get_Y(meas_p, ref_A, ref_B)
    meas_p['x_cor'] = get_X(meas_p, ref_A, ref_B)
    
    meas_p['Sel'], meas_p['err'] = get_Sel(meas_p, ref_C)

# 打印计数结果
print_table()

# 计算结果画图
plt.ioff()
plt.style.use('bmh')

figure = plt.figure('Corrdinate Graph')
sp_handler = figure.add_subplot()
sp_handler.set_xlim((0, test_field_size['x_size']))
sp_handler.set_ylim((0, test_field_size['y_size']))

plot_ref_p(sp_handler, ref_point_list)
plot_meas_p(sp_handler, meas_data_list)

plt.show()    