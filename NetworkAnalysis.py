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

def channelsParser():                                           # Parses and imports raw data that was provided - Code provided as is by Istevan Seres (ELTE)
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

    return channelsData, opens, closes, lifetime

def networkOverTime():                                          # Produces numberOfNodes, numberOfEdges, degrees, avgDegree and maxDegree over time (current time step is 1008 blocks = approx. 1 week) - Code provided as is bz Istevan Seres (ELTE)
    channelsData, opens, closes, lifetime =channelsParser()
    #1008 or 2016 (1 week or 2 weeks measured in blocks)
    LNCreationBlockHeight = 501337
    lastBlock = 576140
    G = nx.MultiGraph()
    seen_nodes = set()
    numberOfNodes= list()
    numberOfEdges = list()
    avgDegree = list()
    maxDegree = list()
    degrees = []
    degreeDistList = []

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
            degrees.append(G.degree())
            degreeDistList.append(nx.degree_histogram(G))
            avgDegree.append(nx.number_of_edges(G)/nx.number_of_nodes(G))
            instance = sorted(G.degree, key=lambda x: x[1], reverse=True)
            maxDegree.append(instance[0][1])

    return  numberOfNodes, numberOfEdges, degrees, degreeDistList, avgDegree, maxDegree

def networkOverTime2():                                         # Produces highdegreee (the Top 10 Nodes with highest Degrees) and networkArray (Network Graph for every week in the Dataset)
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

                G.add_edge(channelsData[j]['from'], channelsData[j]['to'])

        if i in closes:
            closedChannels = closes[i]
            for j in closedChannels:
                if G.has_edge(channelsData[j]['from'],channelsData[j]['to']):
                    G.remove_edge(channelsData[j]['from'],channelsData[j]['to'])
                elif G.has_edge(channelsData[j]['to'],channelsData[j]['from']):
                    G.remove_edge(channelsData[j]['to'],channelsData[j]['from'])
                if G.degree(channelsData[j]['from']) ==0:
                    G.remove_node(channelsData[j]['from'])
                if G.degree(channelsData[j]['to']) ==0:
                    G.remove_node(channelsData[j]['to'])

        if i % 1008 == 0:
            instance = sorted(G.degree, key=lambda x: x[1], reverse=True)
            degreeSorted.append(instance)
            highdegree.append(instance[:10])
            K = deepcopy(G)
            networkArray.append(K)

    return  highdegree, networkArray, degreeSorted

def removeSmallComponents(G, Threshold):                        # Removes all connected_components of graph G with number of nodes smaller than Threshold
    g = deepcopy(G)
    small_components = [component for component in nx.connected_components(g) if len(component) < Threshold]
    for component in small_components:
        g.remove_nodes_from(component)
    return g

def gamma(degreeDistList):                                      # Calculates gamma values over time - weekly
    gamma_list =[]
    for dist in degreeDistList:
        fit = pl.Fit(dist, discrete=True)
        gamma = fit.power_law.alpha
        gamma_list.append(gamma)
    return gamma_list

def clusteringCoefficient (highdegree, networkArray):           # Calculates Clustering Coefficients as well as the Average Clustering Coefficients of the weekly top 10 Nodes in that weekly network
    coefcomplete = [] # List of Lists of Node Addresses and their Clustering Coefficients of he top 10 highest degree nodes over all weekly networks in the dataset

    for i,highestDegreeNodes in enumerate(highdegree):
        weeklyNodes = []

        for node in highestDegreeNodes:
            weeklyNodes.append(nx.clustering(networkArray[i],node[0]))
        coefcomplete.append(weeklyNodes)
        #Average Clustering Coefficient of the top 10 highest degree nodes every 1008 blocks on the btc blockchain
    averageCoef = [np.mean(np.array(coefficients)) for coefficients in coefcomplete]

    return averageCoef, coefcomplete

def precalcBetweennessCentrality (start, step, networkArray):   # Calculates Betweenness Centrality of all Nodes in the current Network weekly and saves these into .json files - start and step denoted in weeks
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        top = nx.betweenness_centrality(g)
        with open(f"betweenness_result{i}.json", 'w') as fp:
            json.dump(top, fp)

def precalcLoadCentrality (start, step, networkArray):          # Calculates Load Centrality of all Nodes in the current Network weekly and saves these into .json files - start and step denoted in weeks
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        top = nx.load_centrality(g)
        with open(f"load_result{i}.json", 'w') as fp:
            json.dump(top, fp)

def betweennessCentrality (start, step, highdegree):            # Extracts the Betweenness Centrality of the top 10 highest degree nodes and averages their values - start and step denoted in weeks
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

    for weekly_values in betweennessCentralityDict.values(): #float import issues for values too small - values approximated as 0
        weekly_values = [value if value is not None else 0 for value in weekly_values]

        averageBetweennessCentrality.append(np.mean(np.array(weekly_values)))

    return betweennessCentralityDict, averageBetweennessCentrality

def loadCentrality (highdegree, networkArray):                  # Extracts the Load Centrality of the top 10 highest degree nodes and averages their values - start and step denoted in weeks
    loadCentralityDict = {}
    averageLoadCentrality = []

    for week in range(start, len(highdegree), step):
        addresses = [node[0] for node in highdegree[week]]
        buffer=[]
        for address in addresses:
            with open(f"load_result{week}.json") as values:
                resultDict = json.load(values)
            buffer.append(resultDict[address])
        loadCentralityDict[week/step] = buffer  #week/step for better indices

    for weekly_values in loadCentralityDict.values(): #float import issues for values too small - values approximated as 0
        weekly_values = [value if value is not None else 0 for value in weekly_values]

        averageLoadCentrality.append(np.mean(np.array(weekly_values)))
    return loadCentralityList, averageLoadCentrality

def averageShortestPathLength (start, step, networkArray):      # Calculates Average Shortest Path Length - start and step denoted in weeks
    averageShortestPathLengthList = []
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        averageShortestPathLengthList.append(nx.average_shortest_path_length(g))

    with open(f"averageShortestPathLengthList{i}.json", 'w') as fp:
        json.dump(averageShortestPathLengthList, fp)

    return averageShortestPathLengthList

def diameter (networkArray, start, step):                       # Calculates Network Diameter - start and step denoted in weeks
    diameterList = []
    for i in range(start, len(networkArray), step):
        g = removeSmallComponents(networkArray[i], 100)
        diameterList.append(nx.diameter(g))

    with open(f"diameter{i}.json", 'w') as fp:
        json.dump(diameterList, fp)

    return diameterList


# Import from Raw Data - Calculation of Coefficients
numberOfNodes, numberOfEdges, degrees, degreeDistList, avgDegree, maxDegree  = networkOverTime()
highdegree, networkArray, degreeSorted = networkOverTime2()
diameterList = diameter(networkArray, 10, 10)
averageShortestPathLengthList = averageShortestPathLength(networkArray, 10, 10)
averageCoef, coefcomplete = clusteringCoefficient(highdegree, networkArray)
precalcBetweennessCentrality(10, 10, networkArray)
betweennessCentralityDict, averageBetweennessCentrality = betweennessCentrality(10, 10, highdegree)


# Import of intermediate results from Json Files - Calculation of Coefficients
diameterList = diameter(networkArray, 10, 10)
averageShortestPathLengthList = averageShortestPathLength(networkArray, 10, 10)
averageCoef, coefcomplete = clusteringCoefficient(highdegree, networkArray)
precalcBetweennessCentrality(10, 10, networkArray)
betweennessCentralityDict, averageBetweennessCentrality = betweennessCentrality(10, 10, highdegree)


# Plotting of coefficients
def plot_number_of_nodes (numberOfNodes):
    x = range(len(numberOfNodes))
    fig = plt.figure()
    plt.plot(numberOfNodes)

    plt.title(r'Number of Nodes over Time')
    plt.xlabel('Time in weeks')
    plt.ylabel(r'$N$')
    fig.savefig('number_of_nodes.png', format='png')
    plt.show()

def plot_number_of_edges (numberOfEdges):
    x = range(len(numberOfEdges))
    fig = plt.figure()
    plt.plot(numberOfEdges)

    plt.title(r'Number of Payment Channels over Time')
    plt.xlabel('Time in weeks')
    plt.ylabel(r'$L$')
    fig.savefig('number_of_edges.png', format='png')
    plt.show()

def plot_average_degree (avgDegree):
    x = range(len(avgDegree))
    fig = plt.figure()
    plt.plot(avgDegree)

    plt.title(r'Average Degree over Time')
    plt.xlabel('Time in weeks')
    plt.ylabel(r'$\langle k \rangle$')
    fig.savefig('average_degree.png', format='png')
    plt.show()

def plot_maximum_degree (maxDegree):
    x = range(len(maxDegree))
    fig = plt.figure()
    plt.plot(maxDegree)

    plt.title(r'Maximum Degree Degree over Time')
    plt.xlabel('Time in weeks')
    plt.ylabel(r'$k_{max}$')
    fig.savefig('max_degree.png', format='png')
    plt.show()

def plot_degree_dist (degreeDistlist, points_in_time):
    multiplot_array = points_in_time
    fig = plt.figure()

    for point in points_in_time:
        x = range(len(degreeDistList[point]))
        y = degreeDistlist[point]
        plt.scatter(x, y, s = 5)

    plt.ylabel('Number of Nodes')
    plt.yscale('log')
    plt.yscale('log')
    plt.ylim(0.7, max(y))

    plt.xlabel('Degree')
    plt.xscale('log')
    plt.xlim(0.7, max(x)+400)

    plt.title('Degree Distribution')
    plt.legend(points_in_time)
    fig.savefig('degree_dist.png', format='png')
    plt.show()

def plot_gamma (gamma_list):
    x = range(len(gamma_list))
    fig = plt.figure()
    plt.plot(gamma_list)

    plt.title(r'Degree Exponent '+r'$\gamma$'+' over time')
    plt.xlabel('Time in weeks')
    plt.ylabel(r'$\gamma$')
    fig.savefig('gamma_over_time.png', format='png')
    plt.show()

def plot_clustering (complete_list_of_cc,average_cc):
    coefcomplete = complete_list_of_cc
    averageCoef = average_cc
    ### Clustering Coefficient Plots
    fig, ax = plt.subplots()
    t = range(len(coefcomplete))
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

def plot_betweenness_centrality (betweennessCentralityDict, averageBetweennessCentrality):
    fig, ax = plt.subplots()
    t = [10, 20, 30, 40, 50, 60, 70]
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

def plot_load_centrality (loadCentralityDict, averageLoadCentrality):
    fig, ax = plt.subplots()
    t = [10, 20, 30, 40, 50, 60, 70]
    y1 = loadCentralityDict
    y2 = averageLoadCentrality

    for te in [1, 2, 3, 4, 5, 6, 7]:
        plt.scatter([te*10]*len(loadCentralityDict[te]), loadCentralityDict[te], s=5, c='b')

    plt.plot(t, y2, c='r')

    ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
    ax.set_ylabel('Load Centrality %c_L$', color='black')

    label1 = mlines.Line2D([],[],color = 'b', label = '$c_b$ for node with 10 highest degrees')
    label2 = mlines.Line2D([],[],color ='r', label = 'Average $c_L$')

    ax.legend(handles = [label1, label2], bbox_to_anchor=(0,1.02,1,0.3), loc='lower left', mode='expand', ncol=2)

    plt.savefig('lc.png')
    plt.show()

def plot_diameter (diameterList):
    fig, ax = plt.subplots()
    t = [10, 20, 30, 40, 50, 60 , 70]
    y = diameterList
    plt.plot(t, y, c='r')
    ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
    ax.set_ylabel('Average Shortest Path Lengths', color='black')
    plt.savefig('averageShortestPathLengths.png')
    plt.show()

def plot_average_shortest_path (averageShortestPathLengthList):
    fig, ax = plt.subplots()
    t = [10, 20, 30, 40, 50, 60 , 70]
    y = averageShortestPathLengthList
    plt.plot(t, y, c='r')
    ax.set_xlabel('Weeks passed since 2017, December 22nd 12 CET')
    ax.set_ylabel('Average Shortest Path Lengths', color='black')
    plt.savefig('averageShortestPathLengths.png')
    plt.show()
