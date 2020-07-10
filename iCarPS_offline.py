# coding=UTF-8
def generate_finally_result_file(inputfile, sampleSeqFile, predictfile, finallyresultfile):
    
    f = open(inputfile)
    for eachline in f:
        if eachline[0] == '>':
            TNote = eachline.strip()
        else:
            TSeq = eachline.strip()
    
    seqfile = []
    f = open(sampleSeqFile)
    for eachline in f:
        if eachline[0] == '>':
            str1 = eachline.strip()
        elif eachline[0] == '\n':
            continue
        else:
            seq = eachline.strip()
            str2 = "%s\t%s"%(str1, seq)
            seqfile.append(str2)
    f.close()

    resultDict = dict()
    resultDict["P_Site"] = ""
    resultDict["N_Site"] = ""
    
    predfile = open(predictfile).readlines()
    for i in range(len(predfile)):
        temp = predfile[i].strip().split()
        PorNP = "P_Site" if int(temp[0]) == 1 else "N_Site"
        resultDict[PorNP] += "%s\n"%(seqfile[i])

    g = open(finallyresultfile, 'w')
    g.write("Carbonylation Site, as follows:\n")
    g.write(resultDict["P_Site"] + "\n")

    g.write("Non-carbonylation Site, as follows:\n")
    g.write(resultDict["N_Site"])
    g.close()
    
    return finallyresultfile


import sys,os
import model.featureExtraction as featureExact
import model.obtainTop_xx_feature as topfeature

# the file path of sorted feature file by using fscore
K_fscore = "./model/K_sort.fscore"
P_fscore = "./model/P_sort.fscore"
R_fscore = "./model/R_sort.fscore"
T_fscore = "./model/T_sort.fscore"

# the parameter files : top feature
K_para = "./model/K.para"
P_para = "./model/P.para"
R_para = "./model/R.para"
T_para = "./model/T.para"


#理化性质文件
PCPFile = "./model/nine_physicochemical_properties_of_amino_acid_Stand.txt" 

predtype = sys.argv[1] #pred_type
inSeqFile = sys.argv[2] #原始输入序列文件

file_path = os.getcwd()#获得当前工作目录
prefix_name = file_path + "/output/" + os.path.basename(inSeqFile).split('.')[0]

csv_file = prefix_name + "_fscore_sorted.csv" #fscore 排好序的特征文件
sampleSeqFile = prefix_name + "_seq.txt" #长序列被截断的文件
predict_result_file = prefix_name + "_predresult.txt" #预测结果
finallyresultfile = prefix_name + "_finalresult.txt"#最终在网页上展示的结果


def run_rf_predict(inSeqFile, predtype, PCPFile, Fscorefile, parafile, csv_file):
    featureExact.main(inSeqFile, predtype, PCPFile, Fscorefile)
    topfeature.main(parafile, csv_file)
   
if __name__ == '__main__':

    if predtype == "K":
       run_rf_predict(inSeqFile, "K", PCPFile, K_fscore, K_para, csv_file)
       os.system("java -jar ./model/prediction.jar ./model/K_RandomForest.model " + prefix_name + "_fscore_sorted_top.arff " +  prefix_name +"_predresult.txt")
       finallyresultfile = generate_finally_result_file(inSeqFile, sampleSeqFile, predict_result_file, finallyresultfile)
    
    elif predtype == "P":
       run_rf_predict(inSeqFile, "P", PCPFile, P_fscore, P_para, csv_file)
       os.system("java -jar ./model/prediction.jar ./model/P_RandomForest.model " + prefix_name + "_fscore_sorted_top.arff " +  prefix_name +"_predresult.txt")
       finallyresultfile = generate_finally_result_file(inSeqFile, sampleSeqFile, predict_result_file, finallyresultfile)
    
    elif predtype == "R":
       run_rf_predict(inSeqFile, "R", PCPFile, R_fscore, R_para, csv_file)
       os.system("java -jar ./model/prediction.jar ./model/R_RandomForest.model " + prefix_name + "_fscore_sorted_top.arff " +  prefix_name +"_predresult.txt")
       finallyresultfile = generate_finally_result_file(inSeqFile, sampleSeqFile, predict_result_file, finallyresultfile)
    
    elif predtype == "T":
       run_rf_predict(inSeqFile, "T", PCPFile, T_fscore, T_para, csv_file)
       os.system("java -jar ./model/prediction.jar ./model/T_RandomForest.model " + prefix_name + "_fscore_sorted_top.arff " +  prefix_name +"_predresult.txt")
       finallyresultfile = generate_finally_result_file(inSeqFile, sampleSeqFile, predict_result_file, finallyresultfile)
    