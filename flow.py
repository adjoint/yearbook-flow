import sys

class Edge:
    def __init__(self, u, v, capacity):
        self.u = u
        self.v = v
        self.capacity = capacity
        self.flow = 0 #not used in the end but why not keep it

#find unique nodes from the list of edges
def identifyNodes(edges):
    nodes = []
    for e in edges:
        if e.u not in nodes:
            nodes.append(e.u)
        if e.v not in nodes:
            nodes.append(e.v)
    return nodes
 
#for each pair of nodes, store their capacity in the list, in a soecific order (11, 12, 13...)
#put 0 if there is no edge
def findCapacityList(edges, nodes):
    n = len(nodes)
    b = [0 for x in range(n*n)]
    for e in edges:
        try:
            b[(e.u-1)*n + e.v-1] = e.capacity
        except IndexError:
            print str(e.v) + " " + str(e.u)
    return b
  
#converting from two nodes to one index to use lists
def getID(u,v,n):
    return (u-1)*n + v - 1
   
#convert a flow problem to slack LP
def initializeAbc(edges,source,sink):
    nodes = identifyNodes(edges)
    #print nodes
    n = len(nodes)
    if n<2:
        return "degenerate case"
    #print n
    capacityConstraints = [0 for x in range(n*n)] #n2
    capacityConstraints2 = findCapacityList(edges, nodes)#n2
    conservationConstraints = [0 for x in range(2*n*n)] #2n2
    sumConstraints = [0 for x in range(2*n-4)] #2*(n-2)
    b = capacityConstraints + capacityConstraints2 + conservationConstraints + sumConstraints
    A = [[0 for x in range(4*n*n + 2*n - 4)] for y in range(2*n*n)]
    #add capacity constraint coefficients
    for node in nodes:
        for node2 in nodes:
            ID = getID(node, node2, n)
            A[n*n+ID][ID] = 1
    #add flow from u to v = - flow from v to u constraint coefficients
    #initially I wanted to do that, but the fact that all flows are nonnegative made it difficult
    #I kept the constaints (with 0s) to not have to worry about my structure, though
    for node in nodes:
        for node2 in nodes:
            ID = getID(node, node2, n)
            ID2 = getID(node2, node, n)
            lst = [0 for x in range(4*n*n + 2*n - 4)]
            if ID==ID2:
                lst[ID] = 0#1
            #print str(node) + " " + str(node2) + " " + str(ID) + " " + str(ID2)
            else:
                lst[ID] = 0#1
                lst[ID2] = 0#1
            A.append(lst)
            lst = [0 for x in range(4*n*n + 2*n - 4)]
            if ID==ID2:
                lst[ID] = 0#-1
            else:
                lst[ID] = 0#-1
                lst[ID2] = 0#-1
            A.append(lst)
    #to account for flow conservation I say the differens between flow in and flow out is zero
    #for each node that's not the source or the sink
    for node in nodes:
        if node!=source and node!=sink:
            lst = [0 for x in range(4*n*n + 2*n - 4)]
            for node2 in nodes:
                if node!=node2:
                    ID = getID(node2, node, n)
                    lst[ID] = 1
            for node2 in nodes:
                if node!=node2:
                    ID = getID(node, node2, n)
                    lst[ID] = -1
            A.append(lst)
            #we need two constraints for each flow conservation constraint
            #because the flow conservation constraint is actually an equality
            lst = [0 for x in range(4*n*n + 2*n - 4)]
            for node2 in nodes:
                if node!=node2:
                    ID = getID(node2, node, n)
                    lst[ID] = -1
            for node2 in nodes:
                if node!=node2:
                    ID = getID(node, node2, n)
                    lst[ID] = 1
            A.append(lst)
    c1 = [1 for x in range(n)]
    c2 = [0 for x in range(4*n*n+n-4)] 
    c = c1 + c2 #the total number of variables we have is 4*n*2*n+n-4 
    #for a in A:
        #print a
    #print "printing b"
    #print b
    #print "printing c"
    #print c
    return (A,b,c,n*n)
   
#We assume the problem is feasible because for maxFlow it is.
#Hence this is not a very generalizable implementation.
#A is size (m+n)*(m+n), and c, b are size m+n, so I need to specify an n -- the number of non-basic variables,
#to know what's going on for sure
def initializeSimplex(A,b,c,n):
    l = len(b)
    noOfNonBasicVar = n
    m = l - noOfNonBasicVar
    N = []
    B = []
    for i in range(noOfNonBasicVar):
        N.append(i+1)
    for i in range(m):
        B.append(noOfNonBasicVar+i+1)
    #print "non-basic variables"
    #print N
    #print "basic variables"
    #print B
    return (N,B,A,b,c,0)
   
def pivot(N,B,A,b,c,v,l,e):
    nA = [[0 for x in range(len(B)+len(N))] for y in range(len(N)+len(B))]
    nb = [0 for x in range(len(B)+len(N))]
    nc = [0 for x in range(len(N)+len(B))]
    nb[e-1] = float(b[l-1])/A[l-1][e-1]
    for j in N:
        if j!= e:
            nA[e-1][j-1] = float(A[l-1][j-1])/float(A[l-1][e-1])
    nA[e-1][l-1] = 1./float(A[l-1][e-1])
    for i in B:
        if i!= l:
            nb[i-1] =  b[i-1] - A[i-1][e-1]*nb[e-1]
            for j in N:
                if j!= e:
                    nA[i-1][j-1] = A[i-1][j-1] - A[i-1][e-1]*nA[e-1][j-1]
            nA[i-1][l-1] = -A[i-1][e-1]*nA[e-1][l-1]
    nv = v + c[e-1]*nb[e-1]
    for j in N:
        if j!= e:
            nc[j-1] = c[j-1] - c[e-1]*nA[e-1][j-1]
    nc[l-1] = -c[e-1]*nA[e-1][l-1]
    nN = []
    nB = []
    for i in N:
        if i!= e:
            nN.append(i)
    nN.append(l)
    for i in B:
        if i!= l:
            nB.append(i)
    nB.append(e)
    #print "PIVOT"
    #for a in nA:
        #print a
    #print "PRINTING b"
    #print nb
    #print "PRINTING c"
    #print nc
    #print "PRINTING v"
    #print nv
    #print "end of pivot"
    return (nN,nB,nA,nb,nc,nv)

def simplex(A,b,c,n):
    (N,B,A,b,c,v) = initializeSimplex(A,b,c,n)        
    someIndexJinNhasCgreaterThanZero = True
    while (someIndexJinNhasCgreaterThanZero == True):
        someIndexJinNhasCgreaterThanZero = False
        for e in N:
            if c[e-1] >0:
                #print "ENTERING TO BASIC: " + str(e)
                minDi = sys.maxint
                minL = 0
                someIndexJinNhasCgreaterThanZero = True
                for i in B:
                    di = sys.maxint
                    if A[i-1][e-1] > 0:
                        di = float(b[i-1])/A[i-1][e-1]
                        if di < minDi:
                            minDi = di
                            minL = i
                    else:
                        di = sys.maxint
                if minDi == sys.maxint:
                    return "unbounded"
                else:
                    #print "LEAVING BASIC, ENTERING TO NONBASIC: " + str(minL)
                    (N,B,A,b,c,v) = pivot(N,B,A,b,c,v,minL,e)
                    break
    x = [0 for i in range(len(N))]#+ len(B))]
    #print "PRINTING B"
    #print B
    #print "PRINTING N"
    #print N
    #print "PRINTING b"
    #print b
    for i in range(1, len(N)+1):#+len(B)):
        if i in B:
            x[i-1] = b[i-1]
        else:
            x[i-1] = 0
    #print x
    objectiveValue = v
    for i in range(len(x)):
        objectiveValue += x[i]*c[i]
    print "objectiveValue = " + str(objectiveValue)
    return x
   
#convert an index of a list to pair of vertices (basically an edge)
def getPairOfVertices(index, n):
    u = index / n
    v = index % n
    u += 1
    v += 1
    return (u,v)
   
def FlowLP(edges,source,sink):
    (A,b,c,n) = initializeAbc(edges,source,sink)
    #print (A,b,c)
    nodes = identifyNodes(edges)
    result = simplex(A,b,c,n)
    nOfNodes = len(nodes)
    capacityList = findCapacityList(edges, nodes)
    flow = result[:len(c)] #flow through initial non-basic variables
    for i in range(len(flow)):
        f = flow[i]
        (u,v) = getPairOfVertices(i, nOfNodes)
        if u!=v and capacityList[i] != 0:
            if int(u)!=int(source) and int(u)!= int(sink) and int(v)!=int(source) and int(v)!=int(sink) and capacityList[i] == f: #we're only interested in edges between people and jobs
                print "edge: (" + str(u) + ", " + str(v) + "), flow: " + str(f) + ", capacity: " + str(capacityList[i]) + ", Assign this job"
            else:
                doNothing = 0
                #print "edge: (" + str(u) + ", " + str(v) + "), flow: " + str(f) + ", capacity: " + str(capacityList[i])
    return result[:len(c)]
   
#edges = [Edge(1,2,4), Edge(1,3,2), Edge(2,4,3), Edge(3,4,1)]
#FlowLP(edges,1,4)
   
#edges = [Edge(1,2,2), Edge(2,3,1)]
#FlowLP(edges,1,3)
  
#edges = [Edge(1,2,3)]
#FlowLP(edges,1,2)
  
#edges = [Edge(1,2,5), Edge(1,3,15), Edge(2,4,5), Edge(2,5,5), Edge(3,4,5), Edge(3,5,5), Edge(4,6,15), Edge(5,6,15)]
#FlowLP(edges,1,6)

edges = [Edge(1,2,2), Edge(1,3,1), Edge(1,4,1), Edge(1,5,1), Edge(1,6,1), Edge(1,7,1), Edge(1,8,1), Edge(1,9,1), Edge(1,10,1), Edge(1,11,1), 
Edge(12,23,1), Edge(13,23,1), Edge(14,23,1), Edge(15,23,1), Edge(16,23,1), Edge(17,23,1), Edge(18,23,1), Edge(19,23,1), Edge(20,23,1), Edge(21,23,1), Edge(22,23,1), 
Edge(2,12,1), Edge(2,14,1), Edge(2,16,1), Edge(2,18,1), Edge(2,21,1), 
Edge(3,12,1), Edge(3,16,1), Edge(3,19,1), Edge(3,22,1), 
Edge(4,12,1), Edge(4,15,1), Edge(4,17,1), Edge(4,18,1), Edge(4,22,1), 
Edge(5,16,1), Edge(6,15,1), Edge(6,17,1), Edge(7,14,1), Edge(8,13,1), 
Edge(9,16,1), Edge(9,19,1), Edge(10,20,1), Edge(11,12,1), Edge(11,17,1), Edge(11,22,1)]

FlowLP(edges,1,23)