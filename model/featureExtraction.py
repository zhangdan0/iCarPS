# coding=UTF-8
#############################################
#1.对输入序列进行截断，根据中心位点为K、P、R、T，分为27长度的序列
#2.提取特征，CC和9PCP
#3.根据排好序的fscore文件，对特征文件进行排序
#############################################
def getSampleSequenceFile(inSeqFile, aaLabel):
    """对输入序列进行截断 
        1个 seq 对应 多个27长度的样本seq
        并判断是否包含非法字符
        return 生成的文件名， 源文件的注释， 源文件的序列
    """
    file_path = os.getcwd()#获得当前工作目录
    prefix_name = file_path + "/output/" + os.path.basename(inSeqFile).split('.')[0]
    outFile =  prefix_name + "_seq.txt"
    
    g = open(outFile, 'w')
    f = open(inSeqFile)
    noteLine = ""
    for eachline in f:
        if eachline[0] == ">":
            noteLine = eachline.strip()           
        else:
            if noteLine != "":
                sequence = eachline.strip()
                sequence = sequence.upper() #转换为大写字母               
                #判断序列是否含有非法字符
                eachline = str(sequence)
                pattern1 =  re.findall(r'B',eachline)
                pattern2 =  re.findall(r'J',eachline)
                pattern3 =  re.findall(r'O',eachline)
                pattern4 =  re.findall(r'Z',eachline)
                pattern5 =  re.findall(r'U',eachline)       
                if pattern1 or pattern2 or pattern3 or pattern4 or pattern5:
                    print('ERROR: The sequence of %s contains illegal characters: B, J, O, Z or U. Please check it!'%(noteLine))
                    continue

                for i in range(len(sequence)):
                    if sequence[i] == aaLabel:
                        if i < 13 and i < len(sequence)-13:
                            seq = 'X'*(13-i) + sequence[:i] + sequence[i: i+13+1]
                        elif i > len(sequence)-13 and i>13:
                            seq = sequence[i-13:i] + sequence[i:] + 'X'*(13-(len(sequence)-i)+1)
                        else:
                            seq = sequence[i-13: i+13+1]
            
                        if len(seq)==27:
                            g.write(noteLine)
                            g.write(" Postion:%d is %s\n"%(i+1, aaLabel))
                            g.write(seq+"\n")

    g.close()

    f.close()
    return outFile, noteLine, sequence

def ObtainCCAndPCPFeature(inputfile,PCPFile):
######################################第一部分：提取CC特征#########################################
######################提取序列坐标,生成各组氨基酸的xyz均值以及所有组的xyz均值########################
    lst = []
    num_seq = 0
    f=open(inputfile,'r')
    for line in f:
        if line.startswith('>'):
            num_seq+=1
        else:
            line = line.upper()
            line = line.strip('\n\r')
            seq_length = len(line)
            F1=F2=F3=F4=F5=F6=F7=F8=F9=F10=F11=F12=F13=F14=F15=F16=F17=F18=F19=F20=0
            for i in range(0,len(line)):
                if line[i] == 'A':
                    F1 += 1
                elif line[i] == 'V':
                    F2 += 1
                elif line[i] == 'L':
                    F3 += 1
                elif line[i] == 'I':
                    F4 += 1
                elif line[i] == 'P':
                    F5 += 1
                elif line[i] == 'F':
                    F6 += 1
                elif line[i] == 'W':
                    F7 += 1
                elif line[i] == 'M':
                    F8 += 1
                elif line[i] == 'G':
                    F9 += 1
                elif line[i] == 'S':
                    F10 += 1
                elif line[i] == 'T':
                    F11 += 1
                elif line[i] == 'C':
                    F12 += 1
                elif line[i] == 'Y':
                    F13 += 1
                elif line[i] == 'N':
                    F14 += 1
                elif line[i] == 'Q':
                    F15 += 1
                elif line[i] == 'K':
                    F16 += 1
                elif line[i] == 'R':
                    F17 += 1
                elif line[i] == 'H':
                    F18 += 1
                elif line[i] == 'D':  
                    F19 += 1
                elif line[i] == 'E':
                    F20 += 1
            Type1Xvalue=F1*(-36.6794)+F2*(-66.6159)+F3*(-94.4033)+F4*(-100.0599)+F5*(-84.6835)+F6*(-101.0966)+F7*(-12.5874)+F8*(-111.3198)
            Type1Yvalue=F1*(58.8302)+F2*(62.1625)+F3*(38.8595)+F4*(20.2361)+F5*(29.1446)+F6*(-79.3826)+F7*(-158.3870)+F8*(-32.9407)
            Type1Zvalue=F1*(-55.9682)+F2*(-73.5564)+F3*(-82.4134)+F4*(-82.4134)+F5*(-72.3001)+F6*(-103.7705)+F7*(-128.2684)+F8*(-93.7201)
            Type2Xvalue=F9*(-24.7561)+F10*(-43.6432)+F11*(-25.2154)+F12*(-52.4159)+F13*(-53.3563)+F14*(-44.2471)+F15*(-65.6344)
            Type2Yvalue=F9*(26.2092)+F10*(25.3163)+F11*(51.3147)+F12*(-25.2563)+F13*(-68.7011)+F14*(-45.4289)+F15*(-24.8605)
            Type2Zvalue=F9*(65.8804)+F10*(92.1974)+F11*(104.4787)+F12*(106.3209)+F13*(158.9550)+F14*(115.8828)+F15*(128.2518)
            Type3Xvalue=F16*(-3.6804)+F17*(-2.1627)+F18*(-0.6273)
            Type3Yvalue=F16*(-1.8878)+F17*(-4.4286)+F18*(4.3459)
            Type3Zvalue=F16*(-146.1415)+F17*(-174.1303)+F18*(-155.1379)
            Type4Xvalue=F19*(-15.3220)+F20*(-16.9336)
            Type4Yvalue=F19*(43.3370)+F20*(-47.8954)
            Type4Zvalue=F19*(-124.9110)+F20*(-138.0496)
            Type1Xmean=Type1Xvalue/8
            Type1Ymean=Type1Yvalue/8
            Type1Zmean=Type1Zvalue/8
            Type2Xmean=Type2Xvalue/7
            Type2Ymean=Type2Yvalue/7
            Type2Zmean=Type2Zvalue/7
            Type3Xmean=Type3Xvalue/3
            Type3Ymean=Type3Yvalue/3
            Type3Zmean=Type3Zvalue/3
            Type4Xmean=Type4Xvalue/2
            Type4Ymean=Type4Yvalue/2
            Type4Zmean=Type4Zvalue/2
            AllXmean=(Type1Xvalue+Type2Xvalue+Type3Xvalue+Type4Xvalue)/seq_length
            AllYmean=(Type1Yvalue+Type2Yvalue+Type3Yvalue+Type4Yvalue)/seq_length
            AllZmean=(Type1Zvalue+Type2Zvalue+Type3Zvalue+Type4Zvalue)/seq_length
        
            lst.append(Type1Xmean)
            lst.append(Type1Ymean)
            lst.append(Type1Zmean)
            lst.append(Type2Xmean)
            lst.append(Type2Ymean)
            lst.append(Type2Zmean)
            lst.append(Type3Xmean)
            lst.append(Type3Ymean)
            lst.append(Type3Zmean)
            lst.append(Type4Xmean)
            lst.append(Type4Ymean)
            lst.append(Type4Zmean)
            lst.append(AllXmean)
            lst.append(AllYmean)
            lst.append(AllZmean) 
    cc = np.array(lst).reshape(num_seq,15)  
    f.close()
###########################提取累积几何中心特征###############################
################第一步，将每条序列中每个氨基酸用空间坐标表示##################
    lst1 = []
    lst2 = []
    lst3 = []
    num_seq = 0
    f=open(inputfile,'r')
    for line in f:
        if line.startswith('>'):
            num_seq+=1
        else:
            line = line.upper()        #将碱基序列转换成大写
            line = line.strip('\n\r')
            for i in range(0,len(line)):
                if line[i] == 'A':
                    sx='-36.6794'
                    sy='58.8302'
                    sz='-55.9682'                    
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'V':
                    sx='-66.6159'
                    sy='62.1625'
                    sz='73.5564'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'L':
                    sx='-94.4003'
                    sy='38.8595'
                    sz='-82.4134'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'I':
                    sx='-100.5090'
                    sy='20.2361'
                    sz='-82.4134'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'P':
                    sx='-84.6835'
                    sy='29.1446'
                    sz='-72.3001'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'F':
                    sx='-101.0996'
                    sy='-79.3826'
                    sz='-103.7705'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)

                elif line[i] == 'W':
                    sx='-12.5874'
                    sy='-158.3870'
                    sz='-128.2684'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'M':
                    sx='-111.3198'
                    sy='-32.9407'
                    sz='-93.7201'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'G':
                    sx='-24.7561'
                    sy='26.2092'
                    sz='65.8804'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'S':
                    sx='-43.6432'
                    sy='25.313163'
                    sz='92.1974'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'T':
                    sx='-25.2154'
                    sy='51.3147'
                    sz='104.3787'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'C':
                    sx='-51.4159'
                    sy='-25.2563'
                    sz='106.3209'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'Y':
                    sx='-53.3563'
                    sy='-68.7011'
                    sz='158.9550'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'N':
                    sx='-44.2471'
                    sy='-45.4280'
                    sz='115.8828'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'Q':
                    sx='-65.6344'
                    sy='-24.8605'
                    sz='128.2518'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'K':
                    sx='-3.6804'
                    sy='-1.8878'
                    sz='-146.1415'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'R':
                    sx='-2.1627'
                    sy='-4.4286'
                    sz='-174.1303'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'H':
                    sx='-0.6273'
                    sy='4.3459'
                    sz='-155.1379'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'D':
                    sx='-15.3220'
                    sy='43.3370'
                    sz='-124.9110'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
                elif line[i] == 'E':
                    sx='-16.9336'
                    sy='-47.8954'
                    sz='-138.0496'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                elif line[i] == 'X':
                    sx='0'
                    sy='0'
                    sz='0'
                    lst1.append(sx)
                    lst2.append(sy)
                    lst3.append(sz)
                    
    Dx = np.array(lst1).reshape(num_seq,27)
    Dy = np.array(lst2).reshape(num_seq,27)
    Dz = np.array(lst3).reshape(num_seq,27)
    f.close()
    
##########第二步，将每条序列中的每个坐标分量的累加求和再取均值，计算出各分量的累加坐标平均值################   
    #x分量的累计平均值
    Dax = []
    for i in range(0,Dx.shape[0]):
        s=[0]*Dx.shape[1]
        s[0]=float(Dx[i][0])
        for j in range(1,Dx.shape[1]):
            s[j]=s[j-1]+float(Dx[i][j])
        ax=sum(s)/Dx.shape[1]
        Dax.append(ax)
    Dax = np.array(Dax)
    Dax = Dax.reshape(Dx.shape[0],1)
    
    #y分量的累计平均值
    Day = []
    for i in range(0,Dy.shape[0]):
        s=[0]*Dy.shape[1]
        s[0]=float(Dy[i][0])
        for j in range(1,Dx.shape[1]):
            s[j]=s[j-1]+float(Dy[i][j])
        ay=sum(s)/Dy.shape[1]
        Day.append(ay)
    Day = np.array(Day)
    Day = Day.reshape(Dy.shape[0],1) 
    
    #z分量的累计平均值
    Daz = []
    for i in range(0,Dz.shape[0]):
        s=[0]*Dz.shape[1]
        s[0]=float(Dz[i][0])
        for j in range(1,Dz.shape[1]):
            s[j]=s[j-1]+float(Dz[i][j])
        az=sum(s)/Dz.shape[1]
        Daz.append(az)
    Daz = np.array(Daz) 
    Daz = Daz.reshape(Dz.shape[0],1) 
#################第三步，将上述坐标分量合并，每条序列都用一个18维向量表示###############
    CCFeaValue = np.hstack((cc,Dax,Day,Daz))
    df_CCFeaValue = pd.DataFrame(CCFeaValue)
    df_CCFeaValue.columns = ['Type1Xmean','Type1Ymean','Type1Zmean',
                         'Type2Xmean','Type2Ymean','Type2Zmean',
                         'Type3Xmean','Type3Ymean','Type3Zmean',
                         'Type4Xmean','Type4Ymean','Type4Zmean',
                         'AllXmean','AllYmean','AllZmean',
                         'Leiji_xmean','Leiji_ymean','Leiji_zmean']
######################################第二部分：提取9_PCPs特征##########################################
    nine = pd.read_csv(PCPFile,sep='\t')
    seq=[]
    f=open(inputfile,'r')
    for line in f:
        line=line.strip('\n')
        if line.startswith('>'):
            pass
        else:
            seq.append(line)
    f.close()
    
    #第一种理化性质
    a1 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a1[i,j]= nine['A'][0]
                
            elif  seq[i][j]== 'C':
                a1[i,j]= nine['C'][0]
                
            elif  seq[i][j]== 'D':
                a1[i,j]= nine['D'][0]

            elif  seq[i][j]== 'E':
                a1[i,j]= nine['E'][0]

            elif  seq[i][j]== 'F':
                a1[i,j]= nine['F'][0]

            elif  seq[i][j]== 'G':
                a1[i,j]= nine['G'][0]
     
            elif  seq[i][j]== 'H':
                a1[i,j]= nine['H'][0]
                
            elif  seq[i][j]== 'I':
                a1[i,j]= nine['I'][0]

            elif  seq[i][j]== 'K':
                a1[i,j]= nine['K'][0]

            elif  seq[i][j]== 'L':
                a1[i,j]= nine['L'][0]

            elif  seq[i][j]== 'M':
                a1[i,j]= nine['M'][0]

            elif  seq[i][j]== 'N':
                a1[i,j]= nine['N'][0]

            elif  seq[i][j]== 'P':
                a1[i,j]= nine['P'][0]

            elif  seq[i][j]== 'Q':
                a1[i,j]= nine['Q'][0]

            elif  seq[i][j]== 'R':
                a1[i,j]= nine['R'][0]

            elif  seq[i][j]== 'S':
                a1[i,j]= nine['S'][0]

            elif  seq[i][j]== 'T':
                a1[i,j]= nine['T'][0]
     
            elif  seq[i][j]== 'V':
                a1[i,j]= nine['V'][0]

            elif  seq[i][j]== 'W':
                a1[i,j]= nine['W'][0]

            elif  seq[i][j]== 'Y':
                a1[i,j]= nine['Y'][0]
                
            elif  seq[i][j]== 'X':
                a1[i,j]= 0

    col_name = []
    df1 = pd.DataFrame(a1)
    for j in range(1,28):
        col_name.append('Hydrophobicity_%s' %j)
    df1.columns = col_name


    #第二种理化性质           
    a2 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a2[i,j]= nine['A'][1]
                
            elif  seq[i][j]== 'C':
                a2[i,j]= nine['C'][1]
                
            elif  seq[i][j]== 'D':
                a2[i,j]= nine['D'][1]

            elif  seq[i][j]== 'E':
                a2[i,j]= nine['E'][1]

            elif  seq[i][j]== 'F':
                a2[i,j]= nine['F'][1]

            elif  seq[i][j]== 'G':
                a2[i,j]= nine['G'][1]
     
            elif  seq[i][j]== 'H':
                a2[i,j]= nine['H'][1]

            elif  seq[i][j]== 'I':
                a2[i,j]= nine['I'][1]

            elif  seq[i][j]== 'K':
                a2[i,j]= nine['K'][1]

            elif  seq[i][j]== 'L':
                a2[i,j]= nine['L'][1]

            elif  seq[i][j]== 'M':
                a2[i,j]= nine['M'][1]

            elif  seq[i][j]== 'N':
                a2[i,j]= nine['N'][1]

            elif  seq[i][j]== 'P':
                a2[i,j]= nine['P'][1]

            elif  seq[i][j]== 'Q':
                a2[i,j]= nine['Q'][1]

            elif  seq[i][j]== 'R':
                a2[i,j]= nine['R'][1]

            elif  seq[i][j]== 'S':
                a2[i,j]= nine['S'][1]

            elif  seq[i][j]== 'T':
                a2[i,j]= nine['T'][1]
     
            elif  seq[i][j]== 'V':
                a2[i,j]= nine['V'][1]

            elif  seq[i][j]== 'W':
                a2[i,j]= nine['W'][1]

            elif  seq[i][j]== 'Y':
                a2[i,j]= nine['Y'][1]
                
            elif  seq[i][j]== 'X':
                a2[i,j]= 0
                
    col_name = []
    df2 = pd.DataFrame(a2)
    for j in range(1,28):
        col_name.append('Hydrophilicity_%s' %j)
    df2.columns = col_name  

    #第三种理化性质           
    a3 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a3[i,j]= nine['A'][2]
                
            elif  seq[i][j]== 'C':
                a3[i,j]= nine['C'][2]
                
            elif  seq[i][j]== 'D':
                a3[i,j]= nine['D'][2]

            elif  seq[i][j]== 'E':
                a3[i,j]= nine['E'][2]

            elif  seq[i][j]== 'F':
                a3[i,j]= nine['F'][2]

            elif  seq[i][j]== 'G':
                a3[i,j]= nine['G'][2]
     
            elif  seq[i][j]== 'H':
                a3[i,j]= nine['H'][2]

            elif  seq[i][j]== 'I':
                a3[i,j]= nine['I'][2]

            elif  seq[i][j]== 'K':
                a3[i,j]= nine['K'][2]

            elif  seq[i][j]== 'L':
                a3[i,j]= nine['L'][2]

            elif  seq[i][j]== 'M':
                a3[i,j]= nine['M'][2]

            elif  seq[i][j]== 'N':
                a3[i,j]= nine['N'][2]

            elif  seq[i][j]== 'P':
                a3[i,j]= nine['P'][2]

            elif  seq[i][j]== 'Q':
                a3[i,j]= nine['Q'][2]

            elif  seq[i][j]== 'R':
                a3[i,j]= nine['R'][2]

            elif  seq[i][j]== 'S':
                a3[i,j]= nine['S'][2]

            elif  seq[i][j]== 'T':
                a3[i,j]= nine['T'][2]
     
            elif  seq[i][j]== 'V':
                a3[i,j]= nine['V'][2]

            elif  seq[i][j]== 'W':
                a3[i,j]= nine['W'][2]

            elif  seq[i][j]== 'Y':
                a3[i,j]= nine['Y'][2]
                
            elif  seq[i][j]== 'X':
                a3[i,j]= 0

    col_name = []
    df3 = pd.DataFrame(a3)
    for j in range(1,28):
        col_name.append('Mass_%s' %j)
    df3.columns = col_name  


    #第四种理化性质           
    a4 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a4[i,j]= nine['A'][3]
                
            elif  seq[i][j]== 'C':
                a4[i,j]= nine['C'][3]
                
            elif  seq[i][j]== 'D':
                a4[i,j]= nine['D'][3]

            elif  seq[i][j]== 'E':
                a4[i,j]= nine['E'][3]

            elif  seq[i][j]== 'F':
                a4[i,j]= nine['F'][3]

            elif  seq[i][j]== 'G':
                a4[i,j]= nine['G'][3]
     
            elif  seq[i][j]== 'H':
                a4[i,j]= nine['H'][3]

            elif  seq[i][j]== 'I':
                a4[i,j]= nine['I'][3]

            elif  seq[i][j]== 'K':
                a4[i,j]= nine['K'][3]

            elif  seq[i][j]== 'L':
                a4[i,j]= nine['L'][3]

            elif  seq[i][j]== 'M':
                a4[i,j]= nine['M'][3]

            elif  seq[i][j]== 'N':
                a4[i,j]= nine['N'][3]

            elif  seq[i][j]== 'P':
                a4[i,j]= nine['P'][3]

            elif  seq[i][j]== 'Q':
                a4[i,j]= nine['Q'][3]

            elif  seq[i][j]== 'R':
                a4[i,j]= nine['R'][3]

            elif  seq[i][j]== 'S':
                a4[i,j]= nine['S'][3]

            elif  seq[i][j]== 'T':
                a4[i,j]= nine['T'][3]
     
            elif  seq[i][j]== 'V':
                a4[i,j]= nine['V'][3]

            elif  seq[i][j]== 'W':
                a4[i,j]= nine['W'][3]

            elif  seq[i][j]== 'Y':
                a4[i,j]= nine['Y'][3]
                
            elif  seq[i][j]== 'X':
                a4[i,j]= 0   
                
    col_name = []
    df4 = pd.DataFrame(a4)
    for j in range(1,28):
        col_name.append('pK1_%s' %j)
    df4.columns = col_name  

    #第五种理化性质           
    a5 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a5[i,j]= nine['A'][4]
                
            elif  seq[i][j]== 'C':
                a5[i,j]= nine['C'][4]
                
            elif  seq[i][j]== 'D':
                a5[i,j]= nine['D'][4]

            elif  seq[i][j]== 'E':
                a5[i,j]= nine['E'][4]

            elif  seq[i][j]== 'F':
                a5[i,j]= nine['F'][4]

            elif  seq[i][j]== 'G':
                a5[i,j]= nine['G'][4]
     
            elif  seq[i][j]== 'H':
                a5[i,j]= nine['H'][4]

            elif  seq[i][j]== 'I':
                a5[i,j]= nine['I'][4]

            elif  seq[i][j]== 'K':
                a5[i,j]= nine['K'][4]

            elif  seq[i][j]== 'L':
                a5[i,j]= nine['L'][4]

            elif  seq[i][j]== 'M':
                a5[i,j]= nine['M'][4]

            elif  seq[i][j]== 'N':
                a5[i,j]= nine['N'][4]

            elif  seq[i][j]== 'P':
                a5[i,j]= nine['P'][4]

            elif  seq[i][j]== 'Q':
                a5[i,j]= nine['Q'][4]

            elif  seq[i][j]== 'R':
                a5[i,j]= nine['R'][4]

            elif  seq[i][j]== 'S':
                a5[i,j]= nine['S'][4]

            elif  seq[i][j]== 'T':
                a5[i,j]= nine['T'][4]
     
            elif  seq[i][j]== 'V':
                a5[i,j]= nine['V'][4]

            elif  seq[i][j]== 'W':
                a5[i,j]= nine['W'][4]

            elif  seq[i][j]== 'Y':
                a5[i,j]= nine['Y'][4]
                
            elif  seq[i][j]== 'X':
                a5[i,j]= 0
            
    col_name = []
    df5 = pd.DataFrame(a5)
    for j in range(1,28):
        col_name.append('pK2_%s' %j)
    df5.columns = col_name        

    #第六种理化性质           
    a6 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a6[i,j]= nine['A'][5]
                
            elif  seq[i][j]== 'C':
                a6[i,j]= nine['C'][5]
                
            elif  seq[i][j]== 'D':
                a6[i,j]= nine['D'][5]

            elif  seq[i][j]== 'E':
                a6[i,j]= nine['E'][5]

            elif  seq[i][j]== 'F':
                a6[i,j]= nine['F'][5]

            elif  seq[i][j]== 'G':
                a6[i,j]= nine['G'][5]
     
            elif  seq[i][j]== 'H':
                a6[i,j]= nine['H'][5]

            elif  seq[i][j]== 'I':
                a6[i,j]= nine['I'][5]

            elif  seq[i][j]== 'K':
                a6[i,j]= nine['K'][5]

            elif  seq[i][j]== 'L':
                a6[i,j]= nine['L'][5]

            elif  seq[i][j]== 'M':
                a6[i,j]= nine['M'][5]

            elif  seq[i][j]== 'N':
                a6[i,j]= nine['N'][5]

            elif  seq[i][j]== 'P':
                a6[i,j]= nine['P'][5]

            elif  seq[i][j]== 'Q':
                a6[i,j]= nine['Q'][5]

            elif  seq[i][j]== 'R':
                a6[i,j]= nine['R'][5]

            elif  seq[i][j]== 'S':
                a6[i,j]= nine['S'][5]

            elif  seq[i][j]== 'T':
                a6[i,j]= nine['T'][5]
     
            elif  seq[i][j]== 'V':
                a6[i,j]= nine['V'][5]

            elif  seq[i][j]== 'W':
                a6[i,j]= nine['W'][5]

            elif  seq[i][j]== 'Y':
                a6[i,j]= nine['Y'][5]
                
            elif  seq[i][j]== 'X':
                a6[i,j]= 0

    col_name = []
    df6 = pd.DataFrame(a6)
    for j in range(1,28):
        col_name.append('pI_%s' %j)
    df6.columns = col_name


    #第七种理化性质           
    a7 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a7[i,j]= nine['A'][6]
                
            elif  seq[i][j]== 'C':
                a7[i,j]= nine['C'][6]
                
            elif  seq[i][j]== 'D':
                a7[i,j]= nine['D'][6]

            elif  seq[i][j]== 'E':
                a7[i,j]= nine['E'][6]

            elif  seq[i][j]== 'F':
                a7[i,j]= nine['F'][6]

            elif  seq[i][j]== 'G':
                a7[i,j]= nine['G'][6]
     
            elif  seq[i][j]== 'H':
                a7[i,j]= nine['H'][6]

            elif  seq[i][j]== 'I':
                a7[i,j]= nine['I'][6]

            elif  seq[i][j]== 'K':
                a7[i,j]= nine['K'][6]

            elif  seq[i][j]== 'L':
                a7[i,j]= nine['L'][6]

            elif  seq[i][j]== 'M':
                a7[i,j]= nine['M'][6]

            elif  seq[i][j]== 'N':
                a7[i,j]= nine['N'][6]

            elif  seq[i][j]== 'P':
                a7[i,j]= nine['P'][6]

            elif  seq[i][j]== 'Q':
                a7[i,j]= nine['Q'][6]

            elif  seq[i][j]== 'R':
                a7[i,j]= nine['R'][6]

            elif  seq[i][j]== 'S':
                a7[i,j]= nine['S'][6]

            elif  seq[i][j]== 'T':
                a7[i,j]= nine['T'][6]
     
            elif  seq[i][j]== 'V':
                a7[i,j]= nine['V'][6]

            elif  seq[i][j]== 'W':
                a7[i,j]= nine['W'][6]

            elif  seq[i][j]== 'Y':
                a7[i,j]= nine['Y'][6]
                
            elif  seq[i][j]== 'X':
                a7[i,j]= 0

    col_name = []
    df7 = pd.DataFrame(a7)
    for j in range(1,28):
        col_name.append('Rigidity_%s' %j)
    df7.columns = col_name


    #第八种理化性质           
    a8 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a8[i,j]= nine['A'][7]
                
            elif  seq[i][j]== 'C':
                a8[i,j]= nine['C'][7]
                
            elif  seq[i][j]== 'D':
                a8[i,j]= nine['D'][7]

            elif  seq[i][j]== 'E':
                a8[i,j]= nine['E'][7]

            elif  seq[i][j]== 'F':
                a8[i,j]= nine['F'][7]

            elif  seq[i][j]== 'G':
                a8[i,j]= nine['G'][7]
     
            elif  seq[i][j]== 'H':
                a8[i,j]= nine['H'][7]

            elif  seq[i][j]== 'I':
                a8[i,j]= nine['I'][7]

            elif  seq[i][j]== 'K':
                a8[i,j]= nine['K'][7]

            elif  seq[i][j]== 'L':
                a8[i,j]= nine['L'][7]

            elif  seq[i][j]== 'M':
                a8[i,j]= nine['M'][7]

            elif  seq[i][j]== 'N':
                a8[i,j]= nine['N'][7]

            elif  seq[i][j]== 'P':
                a8[i,j]= nine['P'][7]

            elif  seq[i][j]== 'Q':
                a8[i,j]= nine['Q'][7]

            elif  seq[i][j]== 'R':
                a8[i,j]= nine['R'][7]

            elif  seq[i][j]== 'S':
                a8[i,j]= nine['S'][7]

            elif  seq[i][j]== 'T':
                a8[i,j]= nine['T'][7]
     
            elif  seq[i][j]== 'V':
                a8[i,j]= nine['V'][7]

            elif  seq[i][j]== 'W':
                a8[i,j]= nine['W'][7]

            elif  seq[i][j]== 'Y':
                a8[i,j]= nine['Y'][7]
                
            elif  seq[i][j]== 'X':
                a8[i,j]= 0

    col_name = []
    df8 = pd.DataFrame(a8)
    for j in range(1,28):
        col_name.append('Flexibility_%s' %j)
    df8.columns = col_name

    #第九种理化性质           
    a9 = np.zeros(shape=(len(seq),27))
    for i in range(len(seq)):
        for j in range(27):
            if seq[i][j]== 'A':
                a9[i,j]= nine['A'][8]
                
            elif  seq[i][j]== 'C':
                a9[i,j]= nine['C'][8]
                
            elif  seq[i][j]== 'D':
                a9[i,j]= nine['D'][8]

            elif  seq[i][j]== 'E':
                a9[i,j]= nine['E'][8]

            elif  seq[i][j]== 'F':
                a9[i,j]= nine['F'][8]

            elif  seq[i][j]== 'G':
                a9[i,j]= nine['G'][8]
     
            elif  seq[i][j]== 'H':
                a9[i,j]= nine['H'][8]

            elif  seq[i][j]== 'I':
                a9[i,j]= nine['I'][8]

            elif  seq[i][j]== 'K':
                a9[i,j]= nine['K'][8]

            elif  seq[i][j]== 'L':
                a9[i,j]= nine['L'][8]

            elif  seq[i][j]== 'M':
                a9[i,j]= nine['M'][8]

            elif  seq[i][j]== 'N':
                a9[i,j]= nine['N'][8]

            elif  seq[i][j]== 'P':
                a9[i,j]= nine['P'][8]

            elif  seq[i][j]== 'Q':
                a9[i,j]= nine['Q'][8]

            elif  seq[i][j]== 'R':
                a9[i,j]= nine['R'][8]

            elif  seq[i][j]== 'S':
                a9[i,j]= nine['S'][8]

            elif  seq[i][j]== 'T':
                a9[i,j]= nine['T'][8]
     
            elif  seq[i][j]== 'V':
                a9[i,j]= nine['V'][8]

            elif  seq[i][j]== 'W':
                a9[i,j]= nine['W'][8]

            elif  seq[i][j]== 'Y':
                a9[i,j]= nine['Y'][8]
                
            elif  seq[i][j]== 'X':
                a9[i,j]= 0
    
    col_name = []
    df9 = pd.DataFrame(a9)
    for j in range(1,28):
        col_name.append('Irreplaceability_%s' %j)
    df9.columns = col_name     
    #按列合并得到所有理化性质值,共27*9=243维
    df_9PCP = pd.concat([df1,df2, df3, df4, df5, df6, df7, df8, df9],axis=1)
    
###################################第三部分：最后合并两种特征值，为最终特征向量（18+243 = 261维）
    df_CCAndPCPFeature = df_CCFeaValue.join(df_9PCP) #按列合并
    
    sampleLabels = []    
    f = open(inputfile,'r')
    for eachline in f:
        if eachline[0] == '>':
            sampleLabels.append('1')
        else:
            sampleSeq = eachline.strip()
    f.close()
    
       
    return df_CCAndPCPFeature, sampleLabels
    
def getFscoreSortFeatureVector(sortedfile, proLabels):
    """通过特征排序的文件，获得排序的特征向量的 index向量
    """
    f = open(sortedfile,'r')
    fscore = []
    for eachline in f:
        tmp = eachline.split()
        if tmp[0] == "Rank":
            continue
        fscore.append(tmp[1])
    f.close()

    indexList = list(
        map(lambda a: proLabels.index(a), fscore)
        )

    return indexList
    
import os
import re
import pandas as pd
import numpy as np
#inSeqFile = sys.argv[1] #原始输入序列文件
#aaLabel = sys.argv[2] #pred_type
#PCPFile = sys.argv[3]#理化性质
#sortedfile = sys.argv[4]#.fscore

def main(inSeqFile,aaLabel,PCPFile,sortedfile):
    file_path = os.getcwd()#获得当前工作目录
    prefix_name = file_path + "/output/" + os.path.basename(inSeqFile).split('.')[0]
    outfeaFile = prefix_name + "_fscore_sorted.csv"
    outFile, noteLine, sequence = getSampleSequenceFile(inSeqFile, aaLabel)
    df_CCAndPCPFeature, sampleLabels = ObtainCCAndPCPFeature(outFile,PCPFile)
    FeaVals = df_CCAndPCPFeature.values
    col_name = df_CCAndPCPFeature.columns.values.tolist()
    indexList = getFscoreSortFeatureVector(sortedfile, col_name)
    # 根据F-score排序后的特征index向量，生成csv格式文件，保存在outfeaFile
    g = open(outfeaFile, 'w')
    for i in range(len(sampleLabels)):
        wLine = sampleLabels[i]
        featureList = list(FeaVals[i])
        for j in range(len(indexList)):
            wLine += ",%.6f"%(featureList[indexList[j]])
        g.write(wLine + "\n")

    g.close()
