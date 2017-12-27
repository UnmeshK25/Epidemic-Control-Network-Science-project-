import matplotlib.pyplot as plt 
import os
from os import listdir
from os.path import isfile, join
from pathlib import Path





between_dict = {}
clustering_dict = {}
hubshutdown_dict = {}
jaccard_dict = {}
random_dict = {}
found_jaccard = False
found_hub = False
found_cluster = False
found_betweenness = False
found_random = False
def read_dir():
    current_dir = ''
    with open("PATH.txt") as f:
        for line in f:
            current_dir = line
    return current_dir

def read_betweenness(number, count):
    global between_dict
    global found_betweenness
    between_dict[count] = list()
    file = "Betweenness/Betweenness_"+number+".csv"
    print("\n between file is = ", file)
    my_file = Path(file)
    if my_file.is_file():
        found_betweenness = True
        with open(file, 'r', encoding='utf-8') as f:
            index = 0
            for line in f:
                if index > 0:
                    currentline = line.split(",")
                    between_dict[count].append(int(currentline[1]))
                index +=1
            
def read_clustering(number, count):
    global clustering_dict
    global found_cluster
    clustering_dict[count] = list()
    file = "clustering/clustering_"+number+".csv"
    print("\n between file is = ", file)
    my_file = Path(file)
    if my_file.is_file():
        found_cluster = True
        with open(file, 'r', encoding='utf-8') as f:
            index = 0
            for line in f:
                if index > 0:
                    currentline = line.split(",")
                    clustering_dict[count].append(int(currentline[1]))
                index +=1
            
def read_hubshutdown(number, count):
    global hubshutdown_dict
    global found_hub
    hubshutdown_dict[count] = list()
    file = "hubshutdown/hubshutdown_"+number+".csv"
    my_file = Path(file)
    print("\n between file is = ", file)
    if my_file.is_file():
        found_hub = True
        with open(file, 'r', encoding='utf-8') as f:
            index = 0
            for line in f:
                if index > 0:
                    currentline = line.split(",")
                    hubshutdown_dict[count].append(int(currentline[1]))
                index +=1
            
            
def read_jaccard(number, count):
    global jaccard_dict
    global found_jaccard
    jaccard_dict[count] = list()
    file = "jaccard/jaccard_"+number+".csv"
    print("\n between file is = ", file)
    my_file = Path(file)
    if my_file.is_file():
        found_jaccard = True
        with open(file, 'r', encoding='utf-8') as f:
            index = 0
            for line in f:
                if index > 0:
                    currentline = line.split(",")
                    jaccard_dict[count].append(int(currentline[1]))
                index +=1
            
def read_random(number, count):
    global random_dict
    global found_random
    random_dict[count] = list()
    file = "random/random_"+number+".csv"
    print("\n between file is = ", file)
    my_file = Path(file)
    if my_file.is_file():
        found_random = True
        with open(file, 'r', encoding='utf-8') as f:
            index = 0
            for line in f:
                if index > 0:
                    currentline = line.split(",")
                    random_dict[count].append(int(currentline[1]))
                index +=1

def main():
    os.chdir(read_dir())
    mypath = "Betweenness/"
    fileb = Path("Betweenness/Betweenness_0000.csv")
    filec = Path("clustering/clustering_0000.csv")
    fileh = Path("hubshutdown/hubshutdown_0000.csv")
    filej = Path("jaccard/jaccard_0000.csv")
    filer = Path("random/random_0000.csv")

    if fileb.is_file():
        mypath = "Betweenness/"
    elif filec.is_file():
        mypath = "clustering/"
    elif filer.is_file():
        mypath = "random/"
    elif fileh.is_file():
        mypath = "hubshutdown/"
    elif filej.is_file():
        mypath = "jaccard/"
   
    print(mypath)
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    print(len(onlyfiles))
    count = 0
    for line in onlyfiles:
        read_betweenness(line[-8:-4], count)
        read_clustering(line[-8:-4], count)
        read_hubshutdown(line[-8:-4], count)
        read_jaccard(line[-8:-4], count)
        read_random(line[-8:-4], count)
        count+=1
        
    count = 0
    for line in onlyfiles:
        percentageCancelled = [0 , 6 , 11 , 16 , 21 , 26 , 31 , 36, 41, 46, 51, 56, 61, 66 , 71, 76, 81 , 86 , 91 , 96 , 101]
        if(found_betweenness):
            plt.plot(percentageCancelled , between_dict[count] , label = "Betweenness",linewidth = 2.1)
        if(found_cluster):
            plt.plot(percentageCancelled , clustering_dict[count] , label = "Clustering",linewidth = 2.1)
        if(found_hub):
            plt.plot(percentageCancelled , hubshutdown_dict[count],  label = "HUB Removal",linewidth = 2.1)
        if(found_jaccard):
            plt.plot(percentageCancelled , jaccard_dict[count] , label = "Jaccard", linewidth = 2.1)
        if(found_random):
            plt.plot(percentageCancelled , random_dict[count] , label = "Random",linewidth = 2.1)
        plt.ylabel('Infected node Count')
        plt.xlabel('Percentage of cancelled flight')
        plt.title('Infected based on Strategey')
        plt.legend()
        
        #plt.show()
        plt.savefig('analyze_'+str(count)+'.png', transparent=True)
        plt.close()
        count +=1
    

    
    
main()