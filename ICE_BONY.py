"""
PROJET PERSONNEL ICE
Auteur : Audrey BONY
Date du rendu : 1 Dec. 2019
"""
## Library
import os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import copy

## Work directory & Data import
path=input("Enter the path of your working directory : ")
os.chdir(path)

## Function

# readFastaFile1
# PARAMETERS : the filename (the file is in fasta format)
# RETURNS : a unique string with all the protein sequences concatenated
# @author: ebecker
def readFastaFile1(filename) :

    # opening the file whose name is filename
    fd = open(filename,'r')
    txt = fd.read()
    fd.close()

    # txt contains all the text of the file.
    # fisrt, I want to seperate the proteins, the symbol that starts a new protein is '>'
    seqs = txt.split('>')[1:]
    s = ""

    for seq in seqs :
        lines = seq.split('\n')[1:]
        for line in lines :
            s = s + line
    return(s)

# lire_word
# PARAMETERS : filename (in txt format)
# RETURNS : a list with all the words
def lire_word(filename) :
    word = open(filename, "r")
    word= word.read()
    word=word.split('\n')[1:]

    # we store the words in capital letters to simplify the writing afterwards
    words = [x.upper() for x in word]
    return words

# trouve_mot
# PARAMETERS : list of words and proteom list
# RETURNS : dictionnary with words and occurance of these & list with unfound words
def trouve_mot(word, proteome) :
    mot_trouve={}
    mot_non_trouve=[]
    for w in word :
        c=0
        substring = w
        n = proteome.count(substring)
        c+= n
        if c != 0 :
            # record word found and occurance
            mot_trouve[w] = c
        else:
            # recording of the word not found
            mot_non_trouve.append(w)
    # we indicate the number of words found
    print(len(word)-len(mot_non_trouve), " words found")
    return(mot_trouve,mot_non_trouve)


# tri_dico
# PARAMETERS : dictionnary with words and occurance of these
# RETURNS : 2 lists containing words sorted alphabetically and occurance
def tri_dico(word_found):
    L=[]
    M=[]
    for clef,valeur in word_found.items():
        # occurance
        L.append(int(valeur))
        # Words
        M.append(clef)
    return(L,M)

# occ_max
# PARAMETERS : dictionnary with words and occurance of these
# RETURNS : a list with the word with the highest occurance and its occurance
def occ_max(word_found):
    maxi = max(word_found.values())
    for word in word_found:
        if word_found[word] == maxi:
            occ_max_word=[word,maxi]
    return occ_max_word

# plus_frequent_10
# PARAMETERS : dictionnary with words and occurance of these
# RETURNS : a list with the 10th more frequent word and their occurance
def plus_frequent_10(word_found):
    plus_freq=[]

    # we must do a deepcopy to avoid modifying the word_found dictionary
    cherche=copy.deepcopy(word_found)
    while len(plus_freq)<10:
        w=occ_max(cherche)
        plus_freq.append(w)
        key=w[0]
        print(w[0], ' appears ',w[1],' times' )
        # we remove the word we just found from the copy of the list of words found
        del cherche[key]

    return(plus_freq)

# plot_bar
# PARAMETERS : a list with the 10th more frequent word and their occurance
# RETURNS : Bar plot of the 10 most frequent words in the human proteome
def plot_bar(lplus_freq):
    height=[]
    l=len(lplus_freq)

    # storage of word occurance
    for i in lplus_freq:
        height.append(i[1])

    # Plot parameters
    BarName = [lplus_freq[0][0], lplus_freq[1][0], lplus_freq[2][0], lplus_freq[3][0], lplus_freq[4][0],lplus_freq[5][0],lplus_freq[6][0],lplus_freq[7][0],lplus_freq[8][0],lplus_freq[9][0]]
    x = [1,2,3,4,5,6,7,8,9,10]
    plt.bar(x, height)
    plt.xticks(x, BarName)
    plt.title("The 10 most frequent words in the human proteome")
    plt.show()

# plot_bar_trouve_pas_trouve
# PARAMETERS : list with found words and a liste with not found words
# RETURNS : Bar plot of words found or not found from 1 to 4 letters
def plot_bar_trouve_pas_trouve(w,word_not_found):
    a=0
    b=0
    c=0
    d=0
    e=0
    a1=0
    b1=0
    c1=0
    d1=0
    e1=0
    for i in w:
        if len(i)==1:
            a+=1
        elif len(i)==2:
            b+=1
        elif len(i)==3:
            c+=1
        elif len(i)==4:
            d+=1
        elif len(i)==5:
            e+=1
    for i in word_not_found:
        if len(i)==1:
            a1+=1
        elif len(i)==2:
            b1+=1
        elif len(i)==3:
            c1+=1
        elif len(i)==4:
            d1+=1
        elif len(i)==5:
            e1+=1

    # Plot parameters
    letter = 5
    found=[a,b,c,d,e]
    not_found=[a1,b1,c1,d1,e1]
    x = [1,2,3,4,5,6,7,8,9,10]

    fig, ax = plt.subplots()
    index = np.arange(letter)
    bar_width = 0.3

    l1 = plt.bar(index,found,bar_width, color='r',  label='Word found')
    l2 = plt.bar(index + bar_width, not_found, bar_width, color='gold', label='Word not found')

    plt.ylabel('Occurance')
    plt.title("Words found or not found from 1 to 5 letters")
    plt.xticks(index + bar_width, ('One-letter word', '2-letter word', '3-letter word', '4-letter word', '5-letter word'))
    plt.legend()
    plt.show()

##Scripte
# Importing Data
print('Importing Data...')
proteome=readFastaFile1('human-proteome.fasta')
word = lire_word("english-words.txt")
print('Data successfully imported')

# Determining which word is present in the proteome and which is not
print("Looking for words in the proteome...")
search=trouve_mot(word, proteome)
word_found=search[0]
word_not_found=search[1]
print('\n')

# Occurrences and words are separated into two lists
lsearch=tri_dico(word_found)
occ=lsearch[0]
w=lsearch[1]

# Word appearing the most
occ_max_word=occ_max(word_found)
print('The most found word is ', occ_max_word[0], ', it appears ',occ_max_word[1],' times' )
print('\n')

# The 10 most frequently appearing words
lplus_freq=plus_frequent_10(word_found)

# Plots
plot_bar(lplus_freq)
plot_bar_trouve_pas_trouve(w,word_not_found)
