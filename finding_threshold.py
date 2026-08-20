def gen_initial(n , conditions = 5):
#Generating 5 different initial conditions
    int_cond = []
    for i in range(conditions):
        x0 = np.random.normal(mean , std , n)
        y0 = np.random.normal(mean , std , n)
        #Add a z0 if the system is 3 dimensional.
        state0 = np.concatenate((x0 , y0))
        int_cond.append(state0)
    return int_cond

def one_simulation(k, state0, args) :
#returns the number of extreme events in a system.
    sol = odeint(system_name , state0 , t, args)
    x_sol = sol[ : , : n]
    x_mean = x_sol.mean(axis = 1)
    actual_x_mean = x_mean[t_threshold : ]
    return num_extreme_events(actual_x_mean)

def has_any_events(parallel , k,int_cond, args):
#Parallelising over the 5 conditions, returns True if any one of the simulations has extreme events. 
    num_ext = parallel(delayed(one_simulation)(k, state0, A, a, B, c, t, n ,adj_sum) for state0 in int_cond)
    return any(i !=0 for i in num_ext)

def find_k(int_cond , t,args, k , lim , delta_k) :
#Finds a coupling threshold by increasing k in amounts of deltak ,returns 0 if k crosses lim. 
    event_counter = 0
    with parallel_config(n_jobs=5 ,backend="loky") :
        with Parallel() as parallel :
            while delta_k > 0.0001 and k <= lim :
                has_events = has_any_events(parallel ,int_cond, A, a, b , c,d , g, I, t, n , k ,adj_sum) 
                if not has_events and event_counter == 0 :
                    k += delta_k
                    print(f"Currently at k = {k}")
                elif has_events and event_counter == 0 :
                #After finding an extreme event, reduces deltak to 0.1 of its value and starts searching again in increments of deltak.
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
