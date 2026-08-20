import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.sparse import csr_array
from scipy.integrate import odeint
from scipy.signal import find_peaks
import networkx as nx
from time import perf_counter
from joblib import Parallel , delayed ,parallel_config
import pandas as pd
a, b , c , d , g , k = 1.0 , 3.2 , 1.0 , 5.0 , 0.9 , 0.1
I = []
i =0
while i < 100 :
    I.append(1.0 + (i / (100-1)) * 0.01)
    i += 1
I = np.array(I)
t = np.linspace(0 , 1 * 10 ** 5 , 1 * 10 **5 + 1)
def hindmarsh_rose(state , t , A , a , b , c , d, g , k , I ,  adj_sum) :
    n = A.shape[0]
    x = state[ : n]
    y = state[n : 2 * n ]
    z = state[2 * n  : ]
    coupling = (k/(n-1)) * ( A.dot(x) - (adj_sum) * (x) )
    dx_dt = y + b * (x **2) - a * (x ** 3) + I + (g) * (z) * (x) + coupling
    dy_dt = c - d * (x **2) - y
    dz_dt = x

    return np.concatenate((dx_dt, dy_dt , dz_dt))

def periodic(list) :
    if len(list) <= 3 :
        return False
    diff = np.diff(np.array(list))
    median_diff = np.median(diff)
    if median_diff == 0 :
        return False
    ratio = diff / median_diff
    return  np.all((ratio < 1.05) & (ratio > 0.95))
def lower_band(peaks , threshold, extreme_event) :
    almost_extreme_events = [i for i in peaks if i > threshold - 0.1 and i < threshold]
    if len(almost_extreme_events) == 0 :
        return False
    elif len(extreme_event) / len(almost_extreme_events) < 0.3333 :
        return True
    else :
        return False

def num_extreme_events(actual_x_mean) :
    peaks_indices, _ = find_peaks(actual_x_mean)
    peaks = actual_x_mean[peaks_indices]
    if len(peaks) == 0 :
        return 0
    threshold = max(15 , np.mean(peaks) + 4 * np.std(peaks))
    true = np.array(peaks) > threshold
    extreme_event = peaks[true]
    extreme_event_times = peaks_indices[true]
    if periodic(extreme_event_times) :
        extreme_event = np.array([])
    if lower_band(peaks , threshold , extreme_event) :
        extreme_event = np.array([])
    if len(extreme_event) / len(peaks) > 0.05 :
        extreme_event = np.array([])
    
    return len(extreme_event)

def gen_initial(n , conditions = 5):
    int_cond = []
    for i in range(conditions):
        x0 = np.random.normal(0.5 , 0.5 , n)
        y0 = np.random.normal(0.2 , 0.5 , n)
        z0 = np.random.normal(0.1 , 0.5 , n)
        state0 = np.concatenate((x0 , y0 ,z0))
        int_cond.append(state0)
    return int_cond

def one_simulation(A, a, b, c, d , g, I , t, n, k, state0, adj_sum) :
    sol = odeint(hindmarsh_rose , state0 , t, args=(A , a , b , c , d , g, I, k, adj_sum) )
    y_sol = (-1) * (sol[ : ,n : 2 * n])
    y_mean = y_sol.mean(axis = 1)
    actual_y_mean = y_mean[20000 : ]
    return num_extreme_events(actual_y_mean)

def has_any_events(parallel ,int_cond, A, a, b, c, d , g, I , t, n , k , adj_sum):
    try :
        num_ext = parallel(delayed(one_simulation)(A, a, b, c, d, g, I, t, n , k, state0,adj_sum) for state0 in int_cond)
        return any(i !=0 for i in num_ext)
    except Exception as e :
        print(f"Timed out at k = {k}")
        return False
def find_k(int_cond , t, A , a , b , c , d , g , I, n , adj_sum , k=0.1 , lim = 12 , delta_k = 0.05) :
    event_counter = 0
    with parallel_config(n_jobs=5 ,backend="loky"):
        with Parallel() as parallel :
            while delta_k > 0.0001 and k <= lim :
                has_events = has_any_events(parallel,int_cond, A, a, b , c,d , g, I, t, n , k ,adj_sum) 
                if not has_events and event_counter == 0 :
                    k += delta_k
                    print(f"Currently at k = {k}")
                elif has_events and event_counter == 0 :
                    print(f"YAY found an extreme event!")
                    event_counter += 1
                    delta_k *= 0.1
                    k -= delta_k
                elif not has_events and event_counter != 0 :
                    k +=  0.1 * delta_k
                elif has_events and event_counter != 0 :
                    delta_k *=0.1
                    k -= delta_k
    if delta_k < 0.0001 :
        return k + delta_k
    elif k > lim :
        print(f"Coupling threshold has exceeded the maximum limit")
        return 0
def final_parallel(t, a , b , c , d ,g ,I, p) :
    print(f"Currently at edge density = {p}")
    G = nx.watts_strogatz_graph(100 , round(p*(99)) , 0.25 )
    n = G.number_of_nodes()
    l = G.number_of_edges()
    edge_density0 = 2 *(l) / (n * (n-1))
    A = nx.to_numpy_array(G)
    adj_sum = A.sum(axis=1)
    if not nx.is_connected(G):
        return None
    algebraic_connectivity = nx.algebraic_connectivity(G)
    int_cond = gen_initial(n , conditions = 5)
    if p == 0.1 :
        coupling_threshold = find_k(int_cond , t, A , a , b , c , d , g , I, n , adj_sum , k =9 ,lim=13 , delta_k=0.05)
    elif p == 0.2 :
        coupling_threshold = find_k(int_cond , t, A , a , b , c , d , g , I, n , adj_sum , k =5 ,lim=10 , delta_k=0.05)
    elif p ==0.3 :
        coupling_threshold = find_k(int_cond , t, A , a , b , c , d , g , I, n , adj_sum , k =2 ,lim=6 , delta_k=0.05)
    elif p >0.3 :
        coupling_threshold = find_k(int_cond , t, A , a , b , c , d , g , I, n , adj_sum , k =0.1 ,lim=4 , delta_k=0.05)
        print(f"coupling thresholds: {coupling_threshold}")
        print(f"edge_density: {edge_density0}")
        print(f"algebraic_connectivites: {algebraic_connectivity}")           
    return [coupling_threshold , edge_density0 ,algebraic_connectivity]
def parallel_for(t, a , b , c , d ,g ,I , edge_densities) :
    list = Parallel(n_jobs=8)(delayed(final_parallel)(t, a , b , c , d ,g ,I, p) for p in edge_densities)
    coupling_thresholds = []
    edge_density = []
    algebraic_connectivites = []
    for i in list:
        coupling_thresholds.append(i[0])
        edge_density.append(i[1])
        algebraic_connectivites.append(i[2])
        df = pd.DataFrame({"coupling threshold" : coupling_thresholds , "edge density" : edge_density , "algebraic connectivity" : algebraic_connectivites})
        df.to_csv("small_world_hrn.csv" ,index= False)
    return df
def main():
    edge_densities = np.linspace(0.1 , 0.8, 8)
    parallel_for(t, a , b , c , d ,g ,I , edge_densities)
    # fig , axs = plt.subplots(2 , 1)
    # axs[0].scatter(edge_densities , coupling_thresholds)
    # axs[0].set_title("kth vs edge density")
    # axs[0].set_xlabel("e")
    # axs[0].set_ylabel("kth")
    # axs[1].scatter(algebraic_connectivities , coupling_thresholds)
    # axs[1].set_title("kth vs algebraic connectivity")
    # axs[1].set_xlabel("lambda2")
    # axs[1].set_ylabel("kth")
    # plt.tight_layout()
    # plt.show()
if __name__ == "__main__" :
    main()
