from ahocorapy.keywordtree import KeywordTree
import textwrap
import random
import pathlib
from termcolor import colored


U="CTAGTTAG"
V="bvbccvCTAnGTTAGvfqvsdqvqCTAGTTAcGvfdCTACGATAGvvfGTTGTTfdvCTAtggAGsfsdfdCTAdddddddddddAGvbcvbcvb"


print("\n","Text : ",V)
print("Motif :",U)


erreur=input("erreur  : ")
#Pi=textwrap.wrap(U, int(erreur))
Pi=[U[i:i+int(erreur)] for i in range(0, len(U), int(erreur))]
print(Pi)
#Aho-Corasick recherch
kwtree = KeywordTree(case_insensitive=True)
for i in range(0,len(Pi)):
    kwtree.add(Pi[i])
kwtree.finalize()
results = kwtree.search_all(V)
#afichage de tout les occurence
Vals=[]
Keyz=[]
for result in results:
    #print(result)
    Vals.append(result[0])
    Keyz.append(result[1])

dictionary = dict(zip(Keyz, Vals))
print(dictionary,"\n")
def AlignS(s,ran):
    l=[]
    ff=""
    for i in range(len(s)):
        ff=ff+"|"
    #l.append(s)
    l.append(colored(s, ran))
    l.append(colored(ff, ran))
    l.append(colored(s, ran))
    re="\n".join(l)
    with open('DNAseq.txt', 'a+') as the_file:
        the_file.writelines(re)
        the_file.writelines("\n")

def AlignSSS(s):
    ff=""
    l=[]
    for i in range(len(s)):
        ff=ff+" "
    l.append(s)
    l.append(ff)
    l.append(ff)
    re= "\n".join(l)
    with open('DNAseq.txt', 'a+') as the_file:
        the_file.writelines(re)
        the_file.writelines("\n")

def aligmentSeq():
    ze=V[0:Keyz[0]]
    AlignSSS(ze)
    foo = ['red', 'green', 'yellow', 'blue', 'magenta','cyan','red', 'green', 'yellow', 'blue', 'magenta','cyan','red', 'green', 'yellow', 'blue', 'magenta','cyan']
    for i in range(0,len(Keyz)-1):
        ran=random.choice(foo)
        gg=V[Keyz[i]:Keyz[i]+len(Vals[i])]
        AlignS(gg,ran)
        dd=V[Keyz[i]+len(Vals[i]):Keyz[i+1]]
        AlignSSS(dd)
    ran = random.choice(foo)
    gg = V[Keyz[len(Keyz)-1]:Keyz[len(Keyz)-1]+len(Vals[len(Keyz)-1])]
    AlignS(gg, ran)
    dd = V[Keyz[len(Keyz)-1]+len(Vals[len(Keyz)-1]):]
    AlignSSS(dd)
aligmentSeq()

def ComineRead(f):
    lines=f.readlines()
    i,j,e=0,1,2
    s,s2,s3="","",""
    while i<=len(lines)-3:
        s=s+lines[i].replace("\n","")
        s2=s2+lines[j].replace("\n","")
        s3=s3+lines[e].replace("\n","")
        i=i+3
        j=j+3
        e=e+3
    return s,s2,s3
readL=open('DNAseq.txt')
r=ComineRead(readL)
print(r[0])
print(r[1])
print(r[2])
readL.close()
file_to_rem = pathlib.Path("DNAseq.txt")
file_to_rem.unlink()