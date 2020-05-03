import networkx as nx
import json
import os.path
import requests
import operator
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import scipy.stats as st
from scipy import stats
from scipy.stats import powerlaw
from matplotlib import pylab
import numpy.random
import powerlaw as pl
from copy import deepcopy

def channelsParser():
    f = open('ln.tsv','r')
    channelsData = dict()
    opens = dict()
    closes = dict()
    counter = 0
    startingBlock = 1000000
    lastBlock = 0
    lifetime = []
    for line in f:
        fields = line.split('\t')
        data = {}
        data['from'] = '0x'+fields[1][3:]
        data['to'] = '0x'+fields[2][3:]
        data['tx'] = '0x'+fields[3][3:]
        data['input'] = fields[4]
        data['amt'] = fields[5]
        data['opened'] = fields[6]
        if int(fields[6])< startingBlock:
            startingBlock = int(fields[6])
        if lastBlock< int(fields[6]):
            lastBlock = int(fields[6])
        data['closed'] = fields[7]
        channelsData[fields[0]]=data
        counter+=1
        if int(fields[6]) in opens:
            opens[int(fields[6])].append(fields[0])
        else:
            opens[int(fields[6])]=[fields[0]]

        if fields[7]=='\\N\n':
            continue
        if int(fields[7]) in closes:
            closes[int(fields[7])].append(fields[0])
        else:
            closes[int(fields[7])]=[fields[0]]
        lifetime.append(int(fields[7]) - int(fields[6]))
    print("First LN channel was created at block height",startingBlock)
    print("We have LN history up to block height", lastBlock)

    return channelsData, opens, closes, lifetime                              #Provided as is by Istevan Seres (ELTE)

def networkOverTime():
    channelsData, opens, closes, lifetime =channelsParser()
    #1008 or 2016 (1 week or 2 weeks measured in blocks)
    LNCreationBlockHeight = 501337
    lastBlock = 576140
    G = nx.MultiGraph()
    seen_nodes = set()
    numberOfNodes= list()
    numberOfEdges = list()
    avgDegree = list()
    effectDiam = list()
    firstConnectedDegrees = []
    for i in range(1500):
        firstConnectedDegrees.append(0)
    for i in range(LNCreationBlockHeight,lastBlock):
        if i in opens:
            createdChannels = opens[i]
            for j in createdChannels:
                if channelsData[j]['from'] not in seen_nodes:
                    G.add_node(channelsData[j]['from'])
                    seen_nodes.add(channelsData[j]['from'])
                if channelsData[j]['to'] not in seen_nodes:
                    G.add_node(channelsData[j]['to'])
                    seen_nodes.add(channelsData[j]['to'])
                if 550000 < i:
                    if type(G.degree(channelsData[j]['from']))==int and type(G.degree(channelsData[j]['to']))==int:
                        firstConnectedDegrees[max(G.degree(channelsData[j]['from']), G.degree(channelsData[j]['to']))] += 1
                    elif type(G.degree(channelsData[j]['from']))==int:
                        firstConnectedDegrees[max(G.degree(channelsData[j]['from']),len(G.degree(channelsData[j]['to'])))] += 1
                    elif type(G.degree(channelsData[j]['to']))==int:
                        firstConnectedDegrees[max(G.degree(channelsData[j]['to']),len(G.degree(channelsData[j]['from'])))] += 1

                G.add_edge(channelsData[j]['from'], channelsData[j]['to'], capacity=channelsData[j]['amt'])
        if i in closes:
            closedChannels = closes[i]
            for j in closedChannels:
                if G.has_edge(channelsData[j]['from'],channelsData[j]['to']):
                    G.remove_edge(channelsData[j]['from'],channelsData[j]['to'])
                else:
                    G.remove_edge(channelsData[j]['to'],channelsData[j]['from'], capacity=channelsData[j]['amt'])
                if G.degree(channelsData[j]['from']) ==0:
                    G.remove_node(channelsData[j]['from'])
                if G.degree(channelsData[j]['to']) ==0:
                    G.remove_node(channelsData[j]['to'])
        if i % 1008 == 0:
            numberOfNodes.append(nx.number_of_nodes(G))
            numberOfEdges.append(nx.number_of_edges(G))
            avgDegree.append(nx.number_of_edges(G)/nx.number_of_nodes(G))

    return  numberOfNodes, numberOfEdges, avgDegree                             #Produces numberOfNodes, numberOfEdges, avgDegree over time (current time step is 1008 blocks = approx. 1 week) - Provided as is bz Istevan Seres (ELTE)

def getFileNames():
    fileNames = []
    for a, b, c in os.walk('../LNdata/lncaptures/lngraph/2019'):
        if (b == []):
            for i in c:
                fileNames.append(str(a) + '/' + i)
    Files = {}
    for i in fileNames:
        f = open(i, "r")
        timestamp = os.path.basename(i)[:-5]
        Files[timestamp] = i
    sortedFiles = sorted(Files.items(), key=lambda t: t[0])
    return sortedFiles

def drawGraph(G):

def averageShortestPathLengths():
    fileNames = getFileNames()
    avgShortestPathLengths=[]
    avgDegree = []
    for i in range(0,len(fileNames)):
        data = readFile(fileNames[i][1])
        G = defineGraph(data)
        #avg degree
        avgDegree.append(G.number_of_edges()/G.order())
        print(i,G.number_of_edges()/G.order())
        #prune smaller components
        toBeDeleted = []
        for j in range(len(list(nx.connected_components(G)))):
            if len(list(nx.connected_components(G))[j]) < 100:
                print(i,j,len(list(nx.connected_components(G))[j]))
                for k in list(nx.connected_components(G))[j]:
                    toBeDeleted.append(k)
        G.remove_nodes_from(toBeDeleted)
        avgShortestPath = nx.algorithms.shortest_paths.generic.average_shortest_path_length(G)
        print(avgShortestPath)
        avgShortestPathLengths.append(avgShortestPath)
        G.clear()

    top = avgShortestPathLengths
    with open('avgShortestPathLengths.json', 'w') as fp:
        json.dump(top, fp)

    fig, ax1 = plt.subplots()
    t = np.arange(0, len(fileNames), 1)
    lns1 = ax1.plot(t, avgShortestPathLengths, 'b-', label='Avg Shortest Paths')
    ax1.set_xlabel('Days passed since 2019, January 30th 12 CET')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('Average Shortest Path Lengths', color='b')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    lns2 = ax2.plot(t, avgDegree, 'r-', label='Average Degree')
    ax2.set_ylabel('Average Out Degree', color='r')
    ax2.tick_params('y', colors='r')

    # added these three lines
    lns = lns1 + lns2# + lns3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='best')

    # Put a nicer background color on the legend.
    # legend.get_frame().set_facecolor('C0')

    fig.tight_layout()
    fig.savefig('averageShortestPathLengths.png')
    plt.show()                  #Provided as is by Istevan Seres (ELTE)
# graph where capacities are tx fees
def defineFeeGraph(data) -> object:
    G = nx.Graph()
    for x in range(len(data)):
        try:
            G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], feebase=data[x]['node1_policy']['fee_base_msat'])
            G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], feerate=data[x]['node1_policy']['fee_rate_milli_msat'])
        except:
            G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], feebase=data[x]['node2_policy']['fee_base_msat'])
            G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], feerate=data[x]['node2_policy']['fee_rate_milli_msat'])
    return G                #Provided as is by Istevan Seres (ELTE)

def defineGraph(data) -> object:
    G = nx.Graph()
    for x in range(len(data)):
        G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], capacity=data[x]['capacity'])
    return G                   #Provided as is by Istevan Seres (ELTE)

def readFile(fileName) -> object:
    with open(fileName) as f:
        data = json.load(f)
    return data['edges']                  #Provided as is by Istevan Seres (ELTE)

###

def removeSmallComponents(G, Threshold):
    g = deepcopy(G)
    small_components = [component for component in nx.connected_components(g) if len(component) < Threshold]
    for component in small_components:
        g.remove_nodes_from(component)
    return g

def networkOverTime2():
    channelsData, opens, closes, lifetime =channelsParser()
    #1008 or 2016 (1 week or 2 weeks measured in blocks)
    LNCreationBlockHeight = 501337
    lastBlock = 576140
    G = nx.Graph()
    networkArray = []
    seen_nodes = set()
    numberOfNodes= list()
    numberOfEdges = list()
    avgDegree = list()
    effectDiam = list()
    firstConnectedDegrees = []
    highdegree = []
    highDegreeNodeAddresses = []
    clusteringCoefficients = []
    degreeSorted = []
    # for i in range(1500):
    #     firstConnectedDegrees.append(0)
    for i in range(LNCreationBlockHeight,lastBlock):
        if i in opens:
            createdChannels = opens[i]
            for j in createdChannels:
                if channelsData[j]['from'] not in seen_nodes:
                    G.add_node(channelsData[j]['from'])
                    seen_nodes.add(channelsData[j]['from'])
                if channelsData[j]['to'] not in seen_nodes:
                    G.add_node(channelsData[j]['to'])
                    seen_nodes.add(channelsData[j]['to'])
                # if 550000 < i:
                #     if type(G.degree(channelsData[j]['from']))==int and type(G.degree(channelsData[j]['to']))==int:
                #         firstConnectedDegrees[max(G.degree(channelsData[j]['from']), G.degree(channelsData[j]['to']))] += 1
                #     elif type(G.degree(channelsData[j]['from']))==int:
                #         firstConnectedDegrees[max(G.degree(channelsData[j]['from']),len(G.degree(channelsData[j]['to'])))] += 1
                #     elif type(G.degree(channelsData[j]['to']))==int:
                #         firstConnectedDegrees[max(G.degree(channelsData[j]['to']),len(G.degree(channelsData[j]['from'])))] += 1

                G.add_edge(channelsData[j]['from'], channelsData[j]['to']) #,capacity=channelsData[j]['amt']

        if i in closes:
            closedChannels = closes[i]
            for j in closedChannels:
                if G.has_edge(channelsData[j]['from'],channelsData[j]['to']):
                    G.remove_edge(channelsData[j]['from'],channelsData[j]['to'])
                elif G.has_edge(channelsData[j]['to'],channelsData[j]['from']):
                    G.remove_edge(channelsData[j]['to'],channelsData[j]['from'])#, capacity=channelsData[j]['amt']
                if G.degree(channelsData[j]['from']) ==0:
                    G.remove_node(channelsData[j]['from'])
                if G.degree(channelsData[j]['to']) ==0:
                    G.remove_node(channelsData[j]['to'])

        if i % 1008 == 0:
            # numberOfNodes.append(nx.number_of_nodes(G))
            # numberOfEdges.append(nx.number_of_edges(G))
            # avgDegree.append(nx.number_of_edges(G)/nx.number_of_nodes(G))

            #List of Lists containing the top 10 highest degree nodes names and their degree for every data point


            instance = sorted(G.degree, key=lambda x: x[1], reverse=True)
            degreeSorted.append(instance)
            highdegree.append(instance[:10])
            K = deepcopy(G)
            networkArray.append(K)

    return  highdegree, networkArray, degreeSorted                             #Produces highdegreee (the Top 10 Nodes with highest Degrees) and networkArray (Network Graph for every week in the Dataset)

def clusteringCoefficient (highdegree, networkArray):
    coefcomplete = [] # List of Lists of Node Addresses and their Clustering Coefficients of he top 10 highest degree nodes over all weekly networks in the dataset

    for i,highestDegreeNodes in enumerate(highdegree):
        weeklyNodes = []

        for node in highestDegreeNodes:
            weeklyNodes.append(nx.clustering(networkArray[i],node[0]))
        coefcomplete.append(weeklyNodes)
        #Average Clustering Coefficient of the top 10 highest degree nodes every 1008 blocks on the btc blockchain
    averageCoef = [np.mean(np.array(coefficients)) for coefficients in coefcomplete]

    return averageCoef, coefcomplete

def precalcBetweennessCentrality (start, step, networkArray):
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        top = nx.betweenness_centrality(g)
        with open(f"betweenness_result{i}.json", 'w') as fp:
            json.dump(top, fp)

def precalcLoadCentrality (range, step, networkArray):
    for i in range(0, range, step):
        g = removeSmallComponents(networkArray[i], 100)
        top = nx.load_centrality(g)
        with open(f"load_result{i}.json", 'w') as fp:
            json.dump(top, fp)

def betweennessCentrality (start, step, highdegree):
    betweennessCentralityDict = {}
    averageBetweennessCentrality = []

    for week in range(start, len(highdegree), step):
        addresses = [node[0] for node in highdegree[week]]
        buffer=[]
        for address in addresses:
            with open(f"betweenness_result{week}.json") as values:
                resultDict = json.load(values)
            buffer.append(resultDict[address])
        betweennessCentralityDict[week/step] = buffer  #week/step for better indices

    for weekly_values in betweennessCentralityDict.values():
        weekly_values = [value if value is not None else 0 for value in weekly_values]

        averageBetweennessCentrality.append(np.mean(np.array(weekly_values)))

    return betweennessCentralityDict, averageBetweennessCentrality

# def loadCentrality (highdegree, networkArray):
    centralityList = {}
    for i in range(74):
        if i%10 == 0:
            addresses = [a_tuple[0] for a_tuple in highdegree[i]]
            buffer=[]
            for address in addresses:
                with open('load_result{}.json'.format(i)) as values:
                    resultDict = json.load(values)
                result = resultDict.get(address)
                buffer.append(result)
            b = i/10
            centralityList[b] = buffer

    # Average Centrality - Replaced None with 0 for now
    averageCentrality = []
    for i in range(len(centralityList)):
        set = centralityList[i]
        for j in range(len(set)):
            entry = set[j]
            if entry == None:
                set[j] = 0
        averageCentrality.append(np.mean(np.array(set)))

    return centralityList, averageCentrality

def averageShortestPathLength (networkArray, start, step):
    averageShortestPathLengthList = []
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        averageShortestPathLengthList.append(nx.average_shortest_path_length(g))

    with open(f"averageShortestPathLengthList{i}.json", 'w') as fp:
        json.dump(averageShortestPathLengthList, fp)

    return averageShortestPathLengthList

def diameter (networkArray, start, step):
    diameterList = []
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        diameterList.append(nx.diameter(g))

    with open(f"diameter{i}.json", 'w') as fp:
        json.dump(diameterList, fp)

    return diameterList

def attackingHighDegrees(G):
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    G.remove_nodes_from(list(nx.connected_components(G))[1])
    percolationThreshold = 0
    PP = []  # probability that a random node belongs to the giant component
    p = []  # percolation threshold
    diameter = []
    for x in range(30):
        highestDegrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
        G.remove_node(next(iter(highestDegrees))[0])
        percolationThreshold+=1
        largest_cc = max(nx.connected_components(G), key=len)
        PP.append(len(largest_cc) / originalGiantComponentSize)
        p.append(percolationThreshold)
        shortpahts=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G.subgraph(largest_cc))
        diameter.append(shortpahts)
        print(percolationThreshold, len(largest_cc), shortpahts)
        if (len(largest_cc) < originalGiantComponentSize / 100):
            print("Percolation Threshold: ", percolationThreshold, percolationThreshold/originalGiantComponentSize)
            break
    print(diameter)
    return (p,PP, diameter)

def attackingHighDegreesNew(G):
    original_cc = max(nx.connected_components(G), key=len)
    subgraph = G.subgraph(list(original_cc))
    original_cc_graph = deepcopy(subgraph)
    originalGiantComponentSize = len(original_cc)
    number_of_removed_nodes = 0
    plot_points = []
    currently_largest_cc_graph = original_cc_graph
    component_list = list(nx.connected_components(currently_largest_cc_graph))
    currently_largest_cc = original_cc

    while len(currently_largest_cc) >= 0.01*originalGiantComponentSize:

        #Next largest component is being picked while the others are saved - New Graph Which Graph should be attacked next is chosen
        currently_largest_cc = component_list[0]
        component_buffer = component_list.pop(0)
        currently_largest_cc_graph = deepcopy(original_cc_graph.subgraph(list(currently_largest_cc))) # A new subgraph is generated every iteration

        # Highest Degree Node is removed
        currently_largest_cc_graph2 = nx.Graph(currently_largest_cc_graph)
        sorted_degrees = sorted(currently_largest_cc_graph.degree(), key=lambda x: x[1], reverse=True)
        highest_degree = sorted_degrees[0]
        currently_largest_cc_graph2.remove_node(highest_degree[0])
        number_of_removed_nodes += 1

        #scatter plot points
        point = (len(currently_largest_cc)/originalGiantComponentSize, number_of_removed_nodes/originalGiantComponentSize)
        plot_points.append(point)

        #All components produced by the removal of the highest degree node in currently_largest_cc_graph are added back to the  component_buffer
        component_buffer = list(component_buffer)
        component_buffer.extend(list(nx.connected_components(currently_largest_cc_graph2)))
        component_list = sorted(component_buffer, key=len, reverse=True)

    perculation_threshold = number_of_removed_nodes/originalGiantComponentSize

    return plot_points, number_of_removed_nodes, perculation_threshold

highdegree, networkArray, degreeSorted = networkOverTime2()

plot_points, number_of_removed_nodes, perculation_threshold = attackingHighDegreesNew(networkArray[30])

#Modes: RANDOM, HIGH DEGREE, BETWEENESS CENTRALITY
#def percolationThreshold (networkArray, start, step)
#    for i in range(start, len(networkArray), step)





testArray = deepcopy(networkArray)

p,PP, diameter = attackingHighDegrees(testArray[30])
len(p)

#graphviz_layout
#networkArray2 = networkArray
#networkArray[1] is networkArray[22]
#looking at the size Distribution of components

diameterList = diameter(networkArray, 10, 10)
averageShortestPathLengthList = averageShortestPathLength(networkArray, 10, 10)

precalcBetweennessCentrality(10, 10, networkArray)
averageCoef, coefcomplete = clusteringCoefficient(highdegree, networkArray)
betweennessCentralityDict, averageBetweennessCentrality = betweennessCentrality(10, 10, highdegree)
betweennessCentralityDict[1]
### Clustering Coefficient Plots
fig, ax = plt.subplots()
#t = np.arange(0, len(centralityList), 1)
t = range(len(coefcomplete))
#lns1 = ax1.plot(t, onlycoef, 'b-', label='Clustering Coefficient')
y1 = coefcomplete
y2 = averageCoef

for te in range(len(coefcomplete)):
    plt.scatter([te]*len(coefcomplete[te]), coefcomplete[te], s=5, c='b')

plt.plot(t, y2, c='r')

ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
ax.set_ylabel('$Clustering Coefficient C$', color='black')

label1 = mlines.Line2D([],[],color = 'b', label = '$C$ for node with 10 highest degrees')
label2 = mlines.Line2D([],[],color ='r', label = 'Average $C$')

ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)

plt.savefig('clustering.png')
plt.show()
### betweennessCentrality plots
fig, ax = plt.subplots()
#t = np.arange(0, len(centralityList), 1)
t = [10, 20, 30, 40, 50, 60, 70]
#lns1 = ax1.plot(t, onlycoef, 'b-', label='Clustering Coefficient')
y1 = betweennessCentralityDict
y2 = averageBetweennessCentrality

for te in [1, 2, 3, 4, 5, 6, 7]:
    plt.scatter([te*10]*len(betweennessCentralityDict[te]), betweennessCentralityDict[te], s=5, c='b')

plt.plot(t, y2, c='r')

ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
ax.set_ylabel('Betweenness Centrality %c_b$', color='black')

label1 = mlines.Line2D([],[],color = 'b', label = '$c_b$ for node with 10 highest degrees')
label2 = mlines.Line2D([],[],color ='r', label = 'Average $c_b$')

ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)

plt.savefig('bc.png')
plt.show()



fig, ax = plt.subplots()
t = range(len(networkArray))
y = componentSizesAllNetworks
for te in range(len(networkArray)):
    plt.scatter([te]*len(componentSizesAllNetworks[te]), componentSizesAllNetworks[te], s=5, c='b')
ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
ax.set_ylabel('Component Sizes', color='black')
plt.savefig('componentSizesAllNetworks.png')
plt.show()


dump = diameterList
with open('diameterlist.json', 'w') as fp:
     json.dump(dump, fp)



fig, ax = plt.subplots()
#t = np.arange(0, len(centralityList), 1)
t = [10, 20, 30, 40, 50, 60 , 70]
#lns1 = ax1.plot(t, onlycoef, 'b-', label='Clustering Coefficient')
y = averageShortestPathLengthList
plt.plot(t, y, c='r')

ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
ax.set_ylabel('Average Shortest Path Lengths', color='black')

# label1 = mlines.Line2D([],[],color = 'b', label = '$c_B$ for node with 10 highest degrees')
# label2 = mlines.Line2D([],[],color ='r', label = 'Average $c_B$')
#
# ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)

plt.savefig('averageShortestPathLengths.png')
plt.show()

fig, ax = plt.subplots()
#t = np.arange(0, len(centralityList), 1)
t = range(len(networkArray))
#lns1 = ax1.plot(t, onlycoef, 'b-', label='Clustering Coefficient')
y = averageShortestPathLengthList
plt.plot(t, y, c='r')

ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
ax.set_ylabel('Average Shortest Path $a$', color='black')

# label1 = mlines.Line2D([],[],color = 'b', label = '$c_B$ for node with 10 highest degrees')
# label2 = mlines.Line2D([],[],color ='r', label = 'Average $c_B$')
#
# ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)

plt.savefig('averageshortestpath.png')
plt.show()


# with open('result73.json', 'r') as file:
#     result73 = json.load(file)
#
#
# type(result73)
# highdegree[73]
# addresses = [a_tuple[0] for a_tuple in highdegree[73]]

### None Investigation # Potentially to small #seems to happen when opening json
# G = networkArray[0]
# with open('result0.json') as values:
#     resultDict = json.load(values)
# print(resultDict.get('0x03fa8ef83983a453095bd9e28238b9755e6705f5bcf6018d8ea89b194525ec06b1'))

### plot

fig, ax = plt.subplots()
#t = np.arange(0, len(centralityList), 1)
t = [0, 10, 20, 30, 40, 50, 60 ,70]
#lns1 = ax1.plot(t, onlycoef, 'b-', label='Clustering Coefficient')
y1 = centralityList
y2 = averageCentrality
#for te, y1e in zip(t, y1):
for te in range(len(centralityList)):
    plt.scatter([te*10]*len(centralityList[te]), centralityList[te], s=1, c='b')

plt.plot(t, y2, c='r')
#lns2 = ax1.plot(t, averageCoef, 'g-', label='Average Clustering Coefficient')
ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
# Make the y-axis label, ticks and tick labels match the line color.
ax.set_ylabel('Load Centrality $c_L$', color='black')
#ax1.tick_params('y', colors='b')
#plt.xticks(range(len(averageCoef)))
#ax2 = ax1.twinx()
#lns3 = ax2.plot(t, averageCoef, 'r-', label='Average Clustering Coefficient')
#ax2.set_ylabel('Average Clustering Coefficient', color='r')
#ax2.tick_params('y', colors='r')
label1 = mlines.Line2D([],[],color = 'b', label = '$c_L$ for node with 10 highest degrees')
label2 = mlines.Line2D([],[],color ='r', label = 'Average $c_L$')
ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)
plt.savefig('lc.png')
plt.show()
