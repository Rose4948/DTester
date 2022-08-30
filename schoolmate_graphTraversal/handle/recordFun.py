#encoding:utf-8
import config
filepath =config.getRecordFundFile()



def recordPathCoverrate(iteration,coverage):#路径覆盖率
    name = "pathCoverRate.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration==0:
        f.writelines("=====================================\n")
    f.writelines(str(coverage)+"\n")
    f.close()

def recordPathCover_un(iteration,uncoverlist):#未覆盖路径表
    name = "pathCover_un.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration==0:
        f.writelines("=====================================\n")
    f.writelines(str(uncoverlist)+"\n")
    f.close()

def recordAllM(iteration,M):#记录迭代中所有的评估矩阵信息：每次迭代都不断记录新的评估矩阵
    name = "M.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration==0:
        f.writelines("=====================================\n")
    f.writelines("第" + str(iteration) + "代 " + "\n")
    for l in M:
        f.write(str(l) + "\n")
    f.close()

def recordTestCaseRunTime(lenth_ind,time):#记录测试数据运行时间
    name = "TestCaseRunTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write(str(lenth_ind)+" "+str(time) + "\n")
    f.close()

def recordInitalCaseTime():#记录测试数据运行时间
    name = "TestCaseRunTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("=====================================\n")
    f.close()
def recordPopTime(iteration,time):#记录测试数据运行时间
    name = "TestPopRunTime.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration==0:
        f.writelines("=====================================\n")
    f.write(str(time) + "\n")
    f.close()
def recordPopGenerationTime(iteration,time):#记录测试数据运行时间
    name = "PopGenerationTime.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration==0:
        f.writelines("=====================================\n")
    f.write(str(time) + "\n")
    f.close()

def recordInitalPop(iteration,pop):#记录初始种群
    name = "pop.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration == 0:
        f.writelines("=====================================\n")
    f.writelines("第" + str(iteration) + "代 " + "\n")
    for op in pop:
        f.write(str(op) + "\n")
    f.close()

def recordPopTestdata(iteration,Testdata,V):#记录种群的变量和数据
    name = "popTestdata.txt"
    file = filepath + name
    f = open(file, "a+")
    if iteration == 0:
        f.writelines("=====================================\n")
    f.writelines("第" + str(iteration) + "代 " + "\n")
    for i in range(len(Testdata)):
        op1=V[i]
        op2=Testdata[i]
        f.write(str(op1) + "\n")
        f.write(str(op2) + "\n")
    f.close()
def recordElitistInfo(iteration,elitist_seq,elitist_data,elitist_M,elitist_V,e_fitness_list,elitist_path):#记录局部运行总时间
    name = "ElitistInfo.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("=====================第" + str(iteration) + "代==================" + "\n")
    for i in range(len(elitist_seq)):
        op = elitist_seq[i]
        f.write(str(op) + "\n")
    for i in range(len(elitist_data)):
        op = elitist_data[i]
        f.write(str(op) + "\n")
    for i in range(len(elitist_M)):
        op = elitist_M[i]
        f.write(str(op) + "\n")
    for i in range(len(elitist_V)):
        op = elitist_V[i]
        f.write(str(op) + "\n")
    f.write(str(e_fitness_list) + "\n")
    f.write(str(elitist_path) + "\n")
    f.close()
def recordInitalElitistInfo():#初始化精英种群
    name = "ElitistInfo.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("************************************\n")
    f.close()
#########废弃###########
def recordLocalUserTime(usertime):#记录局部运行总时间
    name = "LocalUserTime.txt"
    file = filepath + name
    f = open(file, "a+")
    f.writelines("总用时" + str(usertime) + "\n")
    f.write("=====================================\n")
    f.close()

def recordLocalUnCoverPath(ucp):#记录局部未覆盖路径
    name = "UncoverPafterLocal.txt"
    file = filepath + name
    f = open(file, "a+")
    for p in ucp:
        f.writelines( str(p) +"\n")
    f.write("=====================================\n")
    f.close()

def recordCoverPath(pathlist):#记录覆盖路径
    name = "SpathCover.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write(str(pathlist) + "\n")
    f.write("================================================" + "\n")
    f.close()

def recordTestCase(i,path,test_fit,coverage):#记录测试用例的情况
    name = "RunningTestCase.txt"
    file = filepath + name
    f = open(file, "a+")
    f.write("第"+str(i)+"个测试用例"+str(path) + "\n")
    f.write("个体覆盖路径情况："+str(test_fit) + "\n")
    f.write("个体覆盖率：" + str(coverage) + "\n")
    f.write("================================================" + "\n")
    f.close()


if __name__ == '__main__':
    pathindex = 1
    count = 9
    # recordSAcount(pathindex,count)