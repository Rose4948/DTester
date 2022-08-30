# encoding:utf-8
import obtain_efsm_info2
import random
import seq_to_script
import cal_data_fit
import datetime


def obtain_var_from_path(SM, currPathT):
    SM.repeatTranVarDict = {}
    SM.repeatTranFuncDict = {}  ####store repeat transition
    SM.currPathTranVarDict = {}
    pathVarType2 = []
    SM.copyPathInfo()
    for tran in currPathT:  # change the varibale name for the same T
        if currPathT.count(tran) > 1:
            SM.repeatTrans(currPathT)
            break
    allVar=SM.pathInputVar(currPathT)  # identify input variable in events relating to the current path
    SM.pathProProcess(currPathT)  # rewrite identical variables, ---stored in self.pathDefVar
    varname = SM.originalDef
    sp=["aperc","bperc","cperc","dperc","total"]  # teacher
    sp1 = ["bday",  "aday"]  # addressbook
    sp2 = [ "bmonth",  "amonth"]  # addressbook
    sp3 = ["byear","ayear"]# addressbook
    sop_int =["pid","ppa","ppat","telephone"]#openconf
    sop_email=["email","email1","email2"]#openconf
    for num in varname:
        if num =="txtReload" or num in sp or num in sop_int:
            pathVarType2.append('int')
        elif num == "duedate" or num == "assigneddate" or num == "gradedate" or num == "wasdate":
            pathVarType2.append('date')
        elif num in sp2:
            pathVarType2.append('month')
        elif num in sp1:
            pathVarType2.append('day')
        elif num in sp3:
            pathVarType2.append('year')
        elif num in sop_email:
            pathVarType2.append('email')
        else:
            pathVarType2.append('string')
    return pathVarType2,allVar

def GenerateRandomString():
    ''' 在三个范围内，随机选择n个字符组成长度long为n的string型数据  '''
    # basechars = [chr(i) for i in x]  #把ASCII转为具体字符
    str =['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',\
        'A','B','C','D','E','F','G','H','I','J','K','L','M' ,'N','O','P','Q','R','S','T','U','V','W','X','Y','Z' ]
    res = ''
    long = 6
    for i in range(long):
        c = random.choice(str)
        res += c
    return res  # 这样是一个变量的值

def GenerateRandomInt():
    ret = random.choice(range(0,8000))
    return ret
def GenerateRandomdate():
    # 创建日期辅助表
    datestart = '2018-01-01'
    dateend = '2021-12-30'
    # 转为日期格式
    datestart = datetime.datetime.strptime(datestart, '%Y-%m-%d')
    dateend = datetime.datetime.strptime(dateend, '%Y-%m-%d')
    date_list = []
    date_list.append(datestart.strftime('%Y-%m-%d'))
    while datestart < dateend:
        # 日期叠加一天
        datestart += datetime.timedelta(days=+1)
        # 日期转字符串存入列表
        date_list.append(datestart.strftime('%Y-%m-%d'))
    # print(date_list)
    ret = random.choice(date_list)
    return ret

def GenerateRandomMonth():
    month = ["January","February","March","April","May","June","July","August",
             "September","October","November","December",]
    ret = random.choice(month)
    return ret
def GenerateRandomYear():
    ret = random.choice(range(1971,2021))
    return ret

def GenerateRandomEmail():
    res = GenerateRandomString()
    ConEmail =["@qq.com","@136.com","@buct.edu.com"]
    res += random.choice(ConEmail)
    return res

def GenerateRandomDay():
    ret = random.randint(1,31)
    return ret


def Initial_data(vartype):
    solve = []
    for num in vartype:
        if num == 'Boolean':
            data = random.randint(0, 1)
            solve.append(data)
        elif num == 'string':
            str = GenerateRandomString()
            solve.append(str)
        elif num == 'date':
            str = GenerateRandomdate()
            solve.append(str)
        elif num == 'day':
            str = GenerateRandomDay()
            solve.append(str)
        elif num == 'month':
            str = GenerateRandomMonth()
            solve.append(str)
        elif num == 'year':
            str = GenerateRandomYear()
            solve.append(str)
        elif num == 'email':
            str = GenerateRandomEmail()
            solve.append(str)
        else:
            str1 = GenerateRandomInt()
            solve.append(str1)
    return solve
def threeCouple(path,allvar,data):
    return 0,1,1

def testProcee(SM, currPathT):
    # type: (object, object) -> object
    coverage = 0
    pathVarType,allVar = obtain_var_from_path(SM, currPathT)  #获得序列上的变量名和变量类型
    test_fit =[]
    if len(SM.originalDef) > 0:  ##There exist input variables on the path
            vecdata = Initial_data(pathVarType)  #为变量产生数据
            #print "变量上的数据：",vecdata
            #print "变量上的变量：", allVar
            #print "执行当前序列:"
            seq_to_script.runcase(SM.TEvent, currPathT, vecdata) # 执行当前序列
            test_fit,coverage = cal_data_fit.array_spath()  # 获取对各条路径的覆盖情况 test_fit，及完全覆盖路径的覆盖率
            data = vecdata
    else:
        print"The path has no variables"
        data = [0]

    return data,test_fit,coverage,allVar


if __name__ == '__main__':
    # SM = obtain_efsm_info2.obtain_efsm()
    # SM.allPathNum()
    # print "%s has %s states and  %s transitions" % (SM.name, len(SM.stateList), len(SM.transitionList))
    # pathT = ['T3', 'T101', 'T102', 'T104']  #来自模型schoolmate_Vtest.txt
    # populationSize = 10
    # testProcee(SM,pathT,)
    GenerateRandomdate()
