"""
Usage:
    command:
    Epidemic_simulator.py -bjrch [--delay=<days>] [--nsim=<n>] airport.txt routes.txt

Flags:
    -b: betweenness strategy
    -j: jaccard strategy
    -c: clustering coefficent strategy.
    -r: random strategy
    -h: hub removal strategy



Option:
    --delay=<days>  days after exposure to apply edge removal strategy
    --nsim=<n>      number of simulations 
"""

import copy
import getopt
import networkx as nx
import os
import random
import sys
import time


GAMA = 0.0
BETA = 0.0

jaccard_list = list()
betweenness_list = list()
clustering_list = list()

def main():
  
    global GAMA
    global BETA
    DELAY = 0 
    NUM_OF_ITERATION = 10
    BETA=0.5
    GAMA=0.3

    opts, args = getopt.getopt(sys.argv[1:], "brcjh", ["delay=",
                                                         "nsim=",
                                                         "Beta=",
                                                         "Gama="]
                                                            )
    

    if len(args) < 2:
        print(__doc__)
        exit()
      

    airport = args[0]
    flight = args[1]

    strategies = list()

    for options, argm in opts:
        if options == "-b":
            strategies.append("Betweenness")           ### Betweeness
        elif options == "-r":
            strategies.append("random")
        elif options == "-c":
            strategies.append("clustering")
        elif options == "-j":
            strategies.append("jaccard")
        elif options == "-h":
            strategies.append("hubshutdown")
        elif options == "--delay":
            DELAY = int(argm)
        elif options == "--nsim":
            NUM_OF_ITERATION = int(argm)
        elif options == "--Beta":
            BETA = float(argm) 
        elif options == "--Gama":
            GAMA = float(argm)  
            
    
    seed = 100
    random.seed(seed)

    
    network = Load_Network(airport, flight)

    
    degrees = network.degree()
    weights = dict()
    for airport, degree in degrees.items():
        weights[airport] = network.out_degree(airport) +\
                           network.in_degree(airport)
    starts_list = list()
    for i in range(0,NUM_OF_ITERATION):
        start = list()
        while len(start) < 10:
             start_airport = weighted_random(weights)
             if start_airport not in start:
                 start.append(start_airport)
        starts_list.append(start)


    
    currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    path_file=open("PATH.txt","w")
    path_file.write(currenttime)
    os.makedirs(currenttime)
    os.chdir(currenttime)

   
    edgepool = network.edges(data=True)
   


    for strategy in strategies:
    
        print("{0} Mode.".format(strategy) )

        index = 0

        current_method = " "
        Edge_cancel_list = list()
        
        os.makedirs(strategy)   #this creates new folder for strategy and inside this folder, the csv file get stored
        
        iteration = 0
        efforts = [0]
        efforts.extend(range(1,101,5))
        for target in starts_list:
                
            
            output_file = open("{0}/{0}_{1}.csv".format(strategy,
                                                        pad_string(iteration,4)
                                                        ),"w")
            
            output_file.write('"Effort","Total_Infected_Nodes"\n')

            max_number_of_edges =0
            for effort in efforts:
                Removed_Edge = list()
                title = "{0} - {1}%".format(strategy, effort/100)
                results = spreading_controlling_disease(network, Removed_Edge, target,strategy,effort,weights,edgepool,
                                    title=title,DELAY=DELAY)
                total_infected = results["Infected"] + results["Recovered"]
                output_file.write("{0},{1}\n".format(effort/100,total_infected))
                
                if total_infected == 1:
                    for remaining_effort in range(effort+5,101,5):
                        output_file.write("{0},{1}\n".format(remaining_effort/100,
                                                              total_infected))
                    break

            iteration += 1
            output_file.close()



def weighted_random(weights):
    number = random.random() * sum(weights.values())
    for k,v in weights.items():
        if number <= v:
            break
        number -= v
    return k

def pad_string(integer, n):


    string = str(integer)

    while len(string) < n:
        string = "0" + string

    return string



def Load_Network(nodes, edges):


    print("Loading network.")
    G = nx.DiGraph()

    print("\tLoading airports", end="")
    sys.stdout.flush()

    with open(nodes, 'r', encoding='utf-8') as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

            G.add_node(int(entries[0]),
                       country=entries[3],
                       name=entries[1], 
                       lat=entries[6],
                       lon=entries[7])


    print("\t\t\t\t\t[Done]")
    
    print("\tLoading routes",end="")
    sys.stdout.flush()

    edge_count = 0
    error_count = 0
    duplicate_count = 0
    number_of_lines = 1
    with open(edges, 'r', encoding="utf-8") as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            try:
                if G.has_edge(int(entries[3]),int(entries[5])):
                    G[int(entries[3])][int(entries[5])]['Beta'] +=1         ### Beta=Number of flights traveling from source airport to destination airport
                    duplicate_count += 1
                else:
                    if number_of_lines > 1:
                        sourceNode = int(entries[3])
                        endNode = int(entries[5])
                        G.add_edge(sourceNode, endNode )
                        edge_count += 1
                        G.edge[sourceNode][endNode]['Beta'] = 1
            except ValueError:
                # The value doesn't exist
                error_count += 1
                pass
            number_of_lines += 1

    print("\tFinding largest subgraph",end="")
    undirected = G.to_undirected()
    subgraphs = list(nx.connected_component_subgraphs(undirected))
    subgraph_nodes = subgraphs[0].nodes()
    to_remove = list()
    for node in G.nodes():
        if node not in subgraph_nodes:
            to_remove.append(node)
    G.remove_nodes_from(to_remove)
    
    print("\tRemoving isolated vertices",end="")
    in_degree = G.in_degree()
    out_degree = G.out_degree()
    to_remove = [n for n in in_degree if (in_degree[n] + in_degree[n] < 1)] 
    G.remove_nodes_from(to_remove)
    
    print("\tCalculating edge weights",end="")
    G = Edge_weights_calculation(G)
    
    print("\tCalculating clustering coefficents",end="")
    cluster_network = nx.Graph(G)
    lcluster = nx.clustering(cluster_network)
    for i,j in G.edges():
        cluster_sum = lcluster[i] + lcluster[j]
        G[i][j]['cluster'] = cluster_sum
    

    return G


def random_cancel_list(edgepool):
    random_list = random.sample(edgepool, len(edgepool))
    return random_list

def betweenness_cancel_list(G, effort):
    global betweenness_list
    Removed_Edge = list()
    max_number_of_edges = 0
    if(len(betweenness_list) > 1):
        if effort != 0:
            max_number_of_edges = int(len(betweenness_list) * (effort/100))-1
            Removed_Edge = betweenness_list[0:max_number_of_edges] 
    else:
    
        betweennesses = nx.edge_betweenness_centrality(G,
                                                       weight="weight")
        betweenness_list = sorted(betweennesses.keys(), 
                            key=lambda k: betweennesses[k], reverse=True)
        if effort != 0:
            max_number_of_edges = int(len(betweenness_list) * (effort/100))-1
            Removed_Edge = betweenness_list[0:max_number_of_edges]
            
    return Removed_Edge

def clustering_cancel_list( G, edgepool , effort ):
    global clustering_list
    Removed_Edge = list()
    max_number_of_edges = 0
    if(len(clustering_list) > 1):
        if effort != 0:
            max_number_of_edges = int(len(clustering_list) * (effort/100))-1
            Removed_Edge = clustering_list[0:max_number_of_edges]  
    else:
        sorted_cluster = sorted(edgepool, key=lambda k: k[2]['cluster'],
                        reverse=True)
        for cluster_item in sorted_cluster:
            if G[cluster_item[0]][cluster_item[1]]['cluster'] < 2:
                if G[cluster_item[0]][cluster_item[1]]['cluster'] > 0:
                    clustering_list.append((cluster_item[0], cluster_item[1]))
                    
        if effort != 0:
            max_number_of_edges = int(len(clustering_list) * (effort/100))-1
            Removed_Edge = clustering_list[0:max_number_of_edges]
             
    return Removed_Edge

def jaccard_cancel_list(nodes, G , effort):
    print("\n \t Calculating jaccard coefficents",end="")
    

    Removed_Edge = list()
    global jaccard_list
    max_number_of_edges = 0
    if(len(jaccard_list) > 1):
        if effort != 0:
            max_number_of_edges = int(len(jaccard_list) * (effort/100))-1
            Removed_Edge = jaccard_list[0:max_number_of_edges]
    else:
        for node in G:
            successors = G.successors(node)
            for successor in successors:
                jaccard = nx.jaccard_coefficient(nx.Graph(G),[(node,successor)])
                for node1, node2, jacc_coff in jaccard:
                    if jacc_coff < 0.4:
                        jaccard_list.append((node,successor))
                    
        if effort != 0:
            max_number_of_edges = int(len(jaccard_list) * (effort/100))-1          
            Removed_Edge = jaccard_list[0:max_number_of_edges]
        
    return Removed_Edge

                        
    
    
def Edge_weights_calculation(input_network):

    
    G = input_network.copy()
    min = 0
    max = 0
    total = 0
    nde1 = ()
    for node in G.nodes():
        successors = G.successors(node)
        for successor in successors:
            if G[node][successor]['Beta'] > max:   ###Beta= Infection probability
                max= G[node][successor]['Beta']
                nde1=(node,successor)
                print(G[node][successor])
    print(nde1)  

    for node in G.nodes():
        successors = G.successors(node)
        weights = dict()

        total_degree = 0
        for successor in successors:
  
            try:
                total_degree += G.out_degree(successor)
            except TypeError:
                
                pass

        for successor in successors:
            successor_degree = G.out_degree(successor)

            try:
                int(successor_degree)
            except TypeError:
                successor_degree = 0

            if total_degree > 0:
                Infection_probability = successor_degree / \
                                           total_degree
            else:
                Infection_probability = 0

            weights[successor] = Infection_probability
        
        largest_weight = 0
        smallest_weight = 2
        for successor, weight in weights.items():
            if weight > largest_weight:
                largest_weight = weight
            elif weight < smallest_weight:
                smallest_weight = weight

        for successor in successors:
            if largest_weight != smallest_weight:
                relative_weight = (weights[successor] - smallest_weight) /\
                                  (largest_weight - smallest_weight)
            else:
                relative_weight = 0
            G[node][successor]['weight'] = relative_weight
            G[node][successor]['Beta'] = (G[node][successor]['Beta'] )/max
            
       

    return G
def spreading_controlling_disease(input_network,vaccination,starts,current_method,effort,weights,edgepool,file_name = "sir.csv", title="",DELAY=0, RECALCULATE=True):

    global GAMA
    global BETA
    
    print("Simulating Disease.")

    network = input_network.copy()
    infect_network= input_network.copy()

   ## f = open(file_name, "w")
    ##f.write("time, s, e, i, r\n")
    infected_list = list()

    sys.stdout.flush()
    for node in network.nodes():
        infect_network.node[node]["status"]="s"
        infect_network.node[node]["color"] = "green"
        
        network.node[node]["status"] =  "s"
        network.node[node]["color"] = "green"
        network.node[node]["days_of_infection"] = 0
    
    for start in starts:
        infected = start
        infect_network.node[infected]["status"]="i"
        infect_network.node[infected]["color"]  = "red"
        network.node[infected]["status"] = "i"
        network.node[infected]["color"]  = "red"
        infected_list.append(start)
        if isinstance(network,nx.DiGraph):
            in_degree = network.in_degree()[infected] 
            out_degree = network.out_degree()[infected]
            degree = in_degree + out_degree
        else:
            degree = network.degree()[infected]

        print("\t",network.node[infected]["name"],"[",degree,"]")


    if vaccination is not None:
        print("\tVaccinated: ", len(vaccination) )
    else: 
        print("\tVaccinated: None")

    Removed_Edge= list()
    vaccination_nodes = list()
    for Days in range(0,99):
 
        if int(Days) == int(DELAY):
            if current_method == "Betweenness":
                vaccination = betweenness_cancel_list(network,effort)
                print("\n between vacinnation  = ",vaccination)
                network.remove_edges_from(vaccination)
            
            elif current_method == "clustering":
                vaccination = clustering_cancel_list( network, edgepool , effort )
                print("\n clustering vacinnation  = ",vaccination)
                network.remove_edges_from(vaccination)
                
            elif current_method == "random":
                vaccination = random_cancel_list(edgepool)
                print("\n clustering vacinnation  = ",vaccination)
                network.remove_edges_from(vaccination)
                
            elif current_method == "hubshutdown":
                removed_nodes_list = sorted(weights.keys(), key=lambda k: weights[k], reverse=True)
                print("effort = ",effort)
                print("close airport = ",removed_nodes_list)
                if effort != 0:
                    max_number_of_nodes = int(len(removed_nodes_list) * (effort/100))-1
                    vaccination_nodes = removed_nodes_list[0:max_number_of_nodes]
                    print(vaccination_nodes)
                    network.remove_nodes_from(vaccination_nodes)
            elif current_method == "jaccard":
                
                print("in jaccard")
                vaccination = jaccard_cancel_list(infected_list, network,effort)
                print("\n jaccard vacinnation  = ",vaccination)
                network.remove_edges_from(vaccination)
            
            elif vaccination is not None:
                print(DELAY,"on step",Days)
                network.remove_edges_from(vaccination)
                
                if RECALCULATE == True:
                    network = Edge_weights_calculation(network)


        
        Susceptible,Exposed,Infected,Recovered = 0,0,0,0
        print("Susceptible,Exposed,Infected,Recovered")
        
                
        for node in network.nodes():
            status = network.node[node]["status"]
            days_of_infection = network.node[node]["days_of_infection"]
            color = network.node[node]["color"]
            
            if status is "i":
               
                if CheckProbability(GAMA):                                                 ### Gama : Infection to Recovery
                    network.node[node]["status"] = "r"
                    status = "r"
                    network.node[node]["color"] = "purple"
                    infected_list.remove(node)
                    
            if status is "e":
               
                if CheckProbability(BETA):      
                    infect_network.node[node]["status"]="i"
                    infect_network.node[node]["color"]  = "red"                                            ### Beta : Exposed to infection
                    network.node[node]["status"] = "i"
                    status = "i"
                    network.node[node]["color"] = "red"
                    infected_list.append(node)

            if status is "i" and days_of_infection >= 11:
        
                network.node[node]["status"] = "r"
                status = "r"
                network.node[node]["color"] = "purple"
                infected_list.remove(node)
            

            if status is "i":
        
                if days_of_infection > 0:
                    victims = network.successors(node)
                    number_infections = 0
                    for victim in victims:
                        infect_status = network.node[victim]["status"]
                        infect = False  
                                


                        if random.uniform(0,1) <= network[node][victim]['Beta']:
                            infect = True
                            number_infections+=1

                        if infect_status == "s" and infect == True:
                            network.node[victim]["status"] = "e"
                            network.node[victim]["days_of_infection"] = 0
                            network.node[victim]["color"] = "#FF6F00"
                network.node[node]["days_of_infection"] += 1
        
        if vaccination is not None:
            for cancelEdge in vaccination:
                if network.node[cancelEdge[0]]["status"] == "r" and network.node[cancelEdge[1]]["status"] == "r":
                    start= cancelEdge[0]
                    end = cancelEdge[1]
                    network.add_edge(start, end)
                    network[start][end]['Beta'] =1
                    vaccination.remove(cancelEdge)
                    print(network.number_of_edges())
            if RECALCULATE == True:
                network = Edge_weights_calculation(network)         
                
        for node in network.nodes():
            status = network.node[node]["status"]
            days_of_infection = network.node[node]["days_of_infection"]
            color = network.node[node]["color"]

            if status is "s":
                
                Susceptible += 1

            if status is "e":
                Exposed += 1

            if status is "v":
                Susceptible += 1

            elif status is "r":

                Recovered += 1

            elif status is "i":
                
                Infected += 1
        print("{0}, {1}, {2}, {3}, {4}".format(Days, Susceptible, Exposed, Infected, Recovered))

        ##printline = "{0}, {1}, {2}, {3}, {4}".format(Days, Susceptible, Exposed, Infected, Recovered)
        ##f.write(printline + "\n")

       

        if Infected is 0:
            break

    graphMLName = effort * 100
    method=current_method
    nx.write_graphml(network, str(method)+'_'+str(graphMLName)+'_final.graphml') 
    nx.write_graphml(infect_network, str(method)+'_'+str(graphMLName)+'_infected.graphml')       

        
    print("\t----------\n\tS: {0}, I: {1}, R: {2}".format(Susceptible,Infected,Recovered))

    return {"Suscceptable":Susceptible,"Infected":Infected, "Recovered":Recovered}

def CheckProbability(prob):
    #takes in a probability and a potential state that it could change to
    #returns true if the state was changed and false otherwise
    if random.uniform(0,1) < prob:
        return True 
    return False

if __name__ == "__main__":
    main()