import Tkinter
import tkFileDialog
import os

#Browse For File Containing reads
def Import_Reads():
    root = Tkinter.Tk()
    root.withdraw() #use to hide tkinter window
    currdir = os.getcwd()
    file = tkFileDialog.askopenfile(parent=root,mode='rb',title='Choose reads file')
    extension = file.name.split(".")[-1]
    return file.name , extension
################################################################################
def Read_Text(Path):
    with open(Path) as f:
        content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    K_D = []
    Notations = []
    Sequences = []
    content = [Notations.append(x.strip()) for x in content]
    Seq = [Sequences.append(x) for x in Notations if not x.isdigit()]
    Paired_notation = Notations[0].split()
    if len(Paired_notation) >= 2:
        K_D.append(Paired_notation[0])
        K_D.append(Paired_notation[1])
        Sequences = Sequences[1:]
    else:
        new = [K_D.append(x) for x in Notations if x.isdigit()]
    return Sequences , K_D
################################################################################
def Read_Fasta (Path):
    "Read Data From FastaQ File , Return Dict of IDs , Sequence"
    Notations = []
    Notations.append(1)
    Sequences = []
    File = open(Path,'r')
    content = File.readlines()
    content = [Sequences.append(line.strip()) for line in content if not line[0] in '@+' if not line.isdigit() and line.isupper()]
    data = dict(zip(content[0::2], content[1::2]))
    return Sequences,Notations
#######################################################################
def Construct_DEBruijnSingle_Graph(Reads):
    DeBruijn_Graph_Dict = {}
    for read in Reads:
        Sequence = []
        if (Is_in_Dict(DeBruijn_Graph_Dict , read[0:(len(read)- 1)]) == 1):
            Seq = DeBruijn_Graph_Dict[read[0:(len(read)- 1)]]
            Seq.append(read[1:])
        else:
            Sequence.append(read[1:])
            DeBruijn_Graph_Dict[read[0:(len(read)- 1)]] = Sequence
    return DeBruijn_Graph_Dict
#####################################################################

#####################################################################
def Is_in_Dict (Reads_Dict , read):
    Check = 0
    for key , value in Reads_Dict.items():
        if read == key:
            Check = 1
    return Check
######################################################################
######################################################################
def Construct_DEBruijnPaired_Graph(Reads):
    Paired_dict = {}
    for read in Reads:
        Keys = []
        Paired_Keys = []
        Paired_Suffix = []
        Keys = read.split("|")
        Paired_Keys.append(Keys[0][0:len(Keys[0]) - 1])
        Paired_Keys.append(Keys[1][0:len(Keys[0]) - 1])
        if (Is_in_Dict(Paired_dict,tuple(Paired_Keys)) == 1):
            Seq = Paired_dict[tuple(Paired_Keys)]
            new = []
            new.append(Keys[0][1:(len(Keys[0]) + 1)])
            new.append(Keys[1][1:(len(Keys[0]) + 1)])
            List = tuple(new)
            Seq.append(List)
            Paired_dict.update({tuple(Paired_Keys) : Seq})
        else:
            Paired_Suffix.append(Keys[0][1:(len(Keys[0]) + 1)])
            Paired_Suffix.append(Keys[1][1:(len(Keys[0]) + 1)])
            new_tuple = []
            new_tuple.append(tuple(Paired_Suffix))
            Paired_dict.update({tuple(Paired_Keys) : new_tuple})
        #################################
    return Paired_dict
######################################################################
"""
This Function Convert The DeBruijn Graph from Reads into Numericals
Input: 1- Original DeBruijn Graph
Output: 1- Numeric DeBruijn Graph  2- List of all Reads
"""
def Convert_Numeric(DeBruijnGraph):
    Reads_List = []
    Elements = []
    Numeric_DeBruijnGraph = {}
    DeBruijnGraphKeys = DeBruijnGraph.keys()
    for i in range(0, len(DeBruijnGraphKeys)):
        idx2 = []
        if not(DeBruijnGraphKeys[i] in Reads_List):
            Reads_List.append(DeBruijnGraphKeys[i])
        Elements = DeBruijnGraph.get(DeBruijnGraphKeys[i])
        idx1 = Reads_List.index(DeBruijnGraphKeys[i]) + 1
        for j in range(0, len(Elements)):
            if not(Elements[j] in Reads_List):
                Reads_List.append(Elements[j])
            idx2.append(Reads_List.index(Elements[j]) + 1)
        Numeric_DeBruijnGraph[idx1] = idx2
    return Numeric_DeBruijnGraph, Reads_List


"""
This Function Get the Start Key From a Numeric Graph
Input: 1- Numeric Graph  2- List of all Read 
Output: 1- Start Key Value  2- Start Key Index
"""
def Get_DeBruijnGraph_Start(Numeric_DeBruijnGraph, Reads_List):
    StartValue = ''
    StartValue_Index = 0
    Keys = Numeric_DeBruijnGraph.keys()
    Values = Dictionary_Values(Numeric_DeBruijnGraph.values())
    for i in range(0, len(Keys)):
        if not(Keys[i] in Values):
            StartValue_Index = Keys[i]
            break
    StartValue = Reads_List[StartValue_Index - 1]
    return StartValue, StartValue_Index

"""
This Function Converts from List of List of Values of Dictionary into one List
HELPER FUNCTION
"""
def Dictionary_Values(DeBruijnGraphValues):
    Values = []
    for i in range(0, len(DeBruijnGraphValues)):
        List = DeBruijnGraphValues[i]
        for j in range(0, len(List)):
            Values.append(List[j])
    return Values


"""
This Function Get all nodes of the graph
HELPER FUNCTION
"""
def Nodes_Assign(Key, Value_List):
    Nodes_List = []
    for i in range(0, len(Value_List)):
        t = Key, Value_List[i]
        Nodes_List.append(t)
    return Nodes_List



"""
This Function Inserts sub-paths into Eulerian Path
HELPER FUNCTION
"""
def Insert_Path(EulireanPath, SubPath, First_itr_Flag, Knot_Key):
    j = 0
    if First_itr_Flag:
        j = 0
    else:
        indices = [i for i, x in enumerate(EulireanPath) if x == Knot_Key]
        j = max(indices)
    for i in range(0, len(SubPath)):
        EulireanPath.insert((j + 1), SubPath[i])
        j += 1
    return EulireanPath


"""
This Function Assigns indeces to insert in it
HELPER FUNCTION
"""
def Unique_Knot_List(Knot_List):
    seen = set()
    seen_add = seen.add
    return [x for x in Knot_List if not (x in seen or seen_add(x))]


"""
This Function Gets the Eulerian Path
Input: 1- Numeric Debruijn Graph  2- Start Key Index
Output: 1- Eulerian Path
"""
def Get_EulireanPath(Numeric_DeBruijnGraph, StartIndex):
    EulireanPath = []
    SubPath = []
    EulireanPath.append(StartIndex)

    Edges_List = Dictionary_Values(Numeric_DeBruijnGraph.values())
    Edges_List.sort()

    Visited_Edges_List = []
    VisitedNodes = []
    UnVisitedNodes = []
    Knot_Key = []
    Visited_Flag = False
    Unvisited_Flag = False
    First_itr_Flag = True

    while True:
        Visited_Edges_List.sort()
        i = 0
        if (Visited_Edges_List == Edges_List and len(UnVisitedNodes) == 0):
            if len(Knot_Key) == 0:
                val = 0
            else:
                val = Knot_Key[0]
            EulireanPath = Insert_Path(EulireanPath, SubPath, First_itr_Flag, val)
            break
        variable = Numeric_DeBruijnGraph.get(StartIndex)
        Nodes_List = Nodes_Assign(StartIndex, Numeric_DeBruijnGraph.get(StartIndex))
        if len(Nodes_List) > 1:
            for j in range(1, len(Nodes_List)):
                if not(Nodes_List[j] in VisitedNodes) and not(Nodes_List[j] in UnVisitedNodes):
                    UnVisitedNodes.append(Nodes_List[j])
                    t = Nodes_List[0]
                    Knot_Key.append(t[0])
            Knot_Key = Unique_Knot_List(Knot_Key)

        if not(Visited_Flag):
            if not(Nodes_List[i] in VisitedNodes) and not(Unvisited_Flag):
                VisitedNodes.append(Nodes_List[0])
                tuple = Nodes_List[0]
                StartIndex = tuple[1]
                SubPath.append(StartIndex)
                Visited_Edges_List.append(tuple[1])
            else:
                Visited_Flag = True

        if Visited_Flag:
            EulireanPath = Insert_Path(EulireanPath, SubPath, First_itr_Flag, Knot_Key[0])
            if not(First_itr_Flag):
                Knot_Key.remove(Knot_Key[0])
            First_itr_Flag = False
            SubPath = []

        if len(UnVisitedNodes) > 0 and Visited_Flag:
            i = 0
            if not (UnVisitedNodes[i] in VisitedNodes):
                VisitedNodes.append(UnVisitedNodes[i])
                tuple = UnVisitedNodes[i]
                StartIndex = tuple[1]
                SubPath.append(StartIndex)
                Visited_Edges_List.append(tuple[1])
                UnVisitedNodes.remove(UnVisitedNodes[i])
                Visited_Flag = False
            else:
                UnVisitedNodes.remove(UnVisitedNodes[i])

    return EulireanPath

"""
This Function Gets the Genome from Eulerian Path (Works for Single Read)
Input: 1- Eulerian Path, Reads List
Output: 1- Genome string
"""
def Get_Genome_Single_Read(EulireanPath, Reads_List):
    String_Path = ''
    for i in range(0, len(EulireanPath)):
        if i == 0:
            String_Path = Reads_List[EulireanPath[i] - 1]
        else:
            Read = Reads_List[EulireanPath[i] - 1]
            String_Path += Read[len(Read) - 1]
    return String_Path



"""
This Function Gets the Genome from Eulerian Path (Works for Pair Read)
Input: 1- Eulerian Path  2-Reads List  3- Length of Read(k)  4- Gap(d)
Output: 1- Genome string, Prefix String, Suffix String
"""
def Get_Genome_Pair_Read(EulireanPath, Reads_List, k, d):
    Prefix_String = ''
    Suffix_String = ''
    Genome = ''
    for i in range(0, len(EulireanPath)):
        if i == 0:
            t = Reads_List[EulireanPath[i] - 1]
            Prefix_String = t[0]
            Suffix_String = t[1]
        else:
            t = Reads_List[EulireanPath[i] - 1]
            Prefix_Read = t[0]
            Suffix_Read = t[1]
            Prefix_String += Prefix_Read[len(Prefix_Read) - 1]
            Suffix_String += Suffix_Read[len(Suffix_Read) - 1]
    Genome = Overlap_Genome(Prefix_String, Suffix_String, k, d)
    return Genome, Prefix_String, Suffix_String

"""
This Function Overlaps Prefix String on Suffix String
HELPER FUNCTION
"""
def Overlap_Genome(Prefix_String, Suffix_String, k, d):
    idx = k + d
    Genome = Prefix_String[:idx]
    str = Suffix_String[idx-1:]
    Genome += Suffix_String[0:]
    return Genome
##########################################################################
#Main
Path , Extention = Import_Reads()   #if ext == 'txt' ---->read text , else read fasta
if (Extention == 'txt'):
    Seq , List = Read_Text(Path)
elif (Extention == 'fastq'):
    Seq,List = Read_Fasta(Path)
if len(List) < 2:
    Single_dict = Construct_DEBruijnSingle_Graph(Seq)
    Numeric_Graph, Reads_List = Convert_Numeric(Single_dict)
    Start, Start_idx = Get_DeBruijnGraph_Start(Numeric_Graph, Reads_List)
    path = Get_EulireanPath(Numeric_Graph, Start_idx)
    path = Get_Genome_Single_Read(path, Reads_List)
    print path
    print len(path)

elif len(List) >= 2 :
    Paired_dict = Construct_DEBruijnPaired_Graph(Seq)
    graph, Reads = Convert_Numeric(Paired_dict)
    Start, idx = Get_DeBruijnGraph_Start(graph, Reads)
    path = Get_EulireanPath(graph, idx)
    Genome, Prefix, Suffix = Get_Genome_Pair_Read(path, Reads, int(List[0]), int(List[1]))
    print Genome
    print len(Genome)
############################################################################