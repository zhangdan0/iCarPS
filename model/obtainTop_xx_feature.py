# coding=UTF-8
def obtain_top_num(para_file):
    """返回top xx
    """
    tmp = open(para_file).readlines()[0]
    tmp = tmp.strip().split("=")

    return int(tmp[1])
    
def generate_csv_file(csv_file, number, outputfile):
    
    g = open(outputfile,'w')

    f = open(csv_file)
    for eachline in f:
        strings = ''
        temp = eachline.strip().split(",")
        for i in range(number+1):
            strings += ",%s"%temp[i]
        g.write(strings.strip(",") + "\n")
    
    
    f.close()
    g.close()
    return outputfile
    
def productArffFile(csvfile, arfffile, number):
    f=open(csvfile)
    f1=open(arfffile,'w')
    f1.write('@relation Carbonylation'+'\n')
    f1.write('\n')
    
    for x in range(1,int(number)+1):
        substr='@attribute Feature%s numeric'%x
        f1.write(substr+'\n') 
        
    f1.write('@attribute label {1,0}'+'\n'+'\n'+'@data'+'\n')

    for line in f:
        line=line.strip().split(',')
        linenew=','.join(line[1:])+','+line[0]
        f1.write(linenew+'\n')
        
    return arfffile
           
import os
def main(para_file, csv_file):
    featureNumber = obtain_top_num(para_file)
    outputfile = os.path.splitext(csv_file)[0] + "_top.csv"
    outputfile = generate_csv_file(csv_file,featureNumber, outputfile)
    arfffile = os.path.splitext(csv_file)[0] + "_top.arff"
    arfffile = productArffFile(outputfile, arfffile, featureNumber)
    
