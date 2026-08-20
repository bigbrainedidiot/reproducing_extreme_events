#The values of the parameters can be accessed from the paper itself.
def fitz_hugh_nagumo(state , t , A , a , B , c , k, adj_sum) :
    n = A.shape[0]
    x = state[ : n]
    y = state[n : ]
    coupling = (k/(n-1)) * ( A.dot(x) - (adj_sum) * (x) )
    dx_dt = x * (a - x) * (x -1) - y + coupling
    dy_dt = B * x - c * y

    return np.concatenate((dx_dt, dy_dt))

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

def lienard_type_oscillator(state , t , A , a , b , F , g , w, k, adj_sum) :
    n = A.shape[0]
    x = state[ : n]
    y = state[n : ]
    coupling = (k/(n-1)) * ( A.dot(y) - (adj_sum) * (y) )
    dx_dt = y 
    dy_dt = -a * x * y - g * x - b * (x **3) + F * np.sin(w * t) + coupling

    return np.concatenate((dx_dt, dy_dt))
  
def rossler_oscillator(state , t , A , a , b , c, k, adj_sum) :
    n = A.shape[0]
    x = state[ : n]
    y = state[n : 2* n ]
    z = state[2 *n : ]
    coupling_x = (k/(n-1)) * ( A.dot(x) - (adj_sum) * (x) )
    coupling_y = (k/(n-1)) * ( A.dot(y) - (adj_sum) * (y) )
    coupling_z = (k/(n-1)) * ( A.dot(z) - (adj_sum) * (z) )
    dx_dt = -y - z + coupling_x 
    dy_dt = x + a * y + coupling_y
    dz_dt = b + z * (x - c) + coupling_z
    return np.concatenate((dx_dt, dy_dt , dz_dt))
