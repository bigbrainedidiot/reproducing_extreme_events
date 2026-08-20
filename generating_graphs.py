#The function is named the same in all files, but is different depending on the network being generated
#Function for generating random graphs
def final_parallel(t, args, p) :
    print(f"Currently at edge density = {p}") #Real time updates
    G = nx.erdos_renyi_graph(100 , p )
    n = G.number_of_nodes()
    A = nx.to_numpy_array(G)
    adj_sum = A.sum(axis=1)
    if not nx.is_connected(G):
        return None
    algebraic_connectivity = nx.algebraic_connectivity(G)
    int_cond = gen_initial(n , conditions = 5)
    coupling_threshold = find_k(int_cond , t, args, adj_sum , k , lim, delta_k)
  #Printing the current data points, incase the whole run takes a really long time
    print(f"coupling thresholds: {coupling_threshold}")
    print(f"edge_density: {p}")
    print(f"algebraic_connectivites: {algebraic_connectivity}") 
    return [coupling_threshold , p ,algebraic_connectivity]
#Function for generating small world graphs
def final_parallel(t, args ,p) :
    print(f"Currently at edge density = {p}")
    G = nx.watts_strogatz_graph(100 , round((p*(99)) , 0.25 ) 
    n = G.number_of_nodes()
    l = G.number_of_edges()
    edge_density0 = 2 *(l) / (n * (n-1))
    A = nx.to_numpy_array(G)
    adj_sum = A.sum(axis=1)
    if not nx.is_connected(G):
        return None
    algebraic_connectivity = nx.algebraic_connectivity(G)
    int_cond = gen_initial(n , conditions = 5)
    coupling_threshold = find_k(int_cond , t, A, args, n, adj_sum ,k, lim , delta_k)
    print(f"coupling thresholds: {coupling_threshold}")
    print(f"edge_density: {edge_density0}")
    print(f"algebraic_connectivites: {algebraic_connectivity}")
#Function for generating scale free graphs
def final_parallel(t, args, p) :
    print(f"Currently at edge density = {p}")
    G = nx.barabasi_albert_graph(100 , round((p*(99))/2) )
    n = G.number_of_nodes()
    l = G.number_of_edges()
    edge_density0 = 2 *(l) / (n * (n-1))
    A = nx.to_numpy_array(G)
    adj_sum = A.sum(axis=1)
    if not nx.is_connected(G):
        return None
    algebraic_connectivity = nx.algebraic_connectivity(G)
    int_cond = gen_initial(n , conditions = 5)
    coupling_threshold = find_k(int_cond , t, A, args, n , adj_sum , k , lim  , delta_k)
    print(f"coupling threshold: {coupling_threshold}")
    print(f"edge_density: {edge_density0}")
    print(f"algebraic_connectivity: {algebraic_connectivity}")  
    return [coupling_threshold , edge_density0 ,algebraic_connectivity]

#Final run
def parallel_for(t, args, edge_densities) :
#Parallelises over 8 edge densities, returns a csv file with all the data after the run
    list = Parallel(n_jobs=8)(delayed(final_parallel)(t, args, p) for p in edge_densities)
    coupling_thresholds = []
    edge_density = []
    algebraic_connectivites = []
    for i in list:
        coupling_thresholds.append(i[0])
        edge_density.append(i[1])
        algebraic_connectivites.append(i[2])
        df = pd.DataFrame({"coupling threshold" : coupling_thresholds , "edge density" : edge_density , "algebraic connectivity" : algebraic_connectivites})
        df.to_csv("file_name.csv" ,index= False)
    return df
def main():
#Calling the main function
    edge_densities = np.linspace(0.1 , 0.8, 8)
    parallel_for(t, args, edge_densities)
