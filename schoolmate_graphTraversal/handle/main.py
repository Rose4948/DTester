# encoding:utf-8
from __future__ import division
import sys

sys.path.append("..")
sys.path.append("../..")
sys.path.append("../module")
from datetime import datetime
import obtain_efsm_info2  # 模型处理器
import generate_seq_ga_check  # 序列生成器-----产生覆盖敏感路径的序列
import execute
import local_search_alg
import recordFun
import config
import seq_to_script

def main():
    LunHui = 0
    fileName_F = config.getRecordFundFile()
    for L in ["ElitistInfo.txt", "popTestdata.txt", "M.txt", "pathCoverRate.txt","pathCover_un.txt", "pop.txt","PopGenerationTime.txt","TestCaseRunTime.txt", "TestPopRunTime.txt","SpathCover.txt"]:
        fileName = fileName_F + L
        f = open(fileName, "r+")
        f.truncate()
        f.close()
    while (LunHui < 3):
        SM = obtain_efsm_info2.obtain_efsm()
        popsize, pc, pm, Max = config.getPopParameter()
        generate_seq_ga_check.mutantion3_Couple_Pre(SM)  #########################new couple
        print "again time：", LunHui
        startTime = datetime.now()
        pop_time0 = datetime.now()
        #pop = generate_seq_ga_check.initialpop_feasible(popsize, SM)
        #control initial pop keeping consistent or use generate_seq_ga_check.initialpop_feasible to randomly generate pop
        pop=[['T96', 'T97', 'T98', 'T109', 'T110', 'T110'],
        ['T125', 'T54', 'T55', 'T55'],
        ['T125', 'T56', 'T94', 'T95', 'T70', 'T71', 'T68', 'T69', 'T72', 'T73', 'T74', 'T68'],
        ['T125', 'T38', 'T41', 'T42', 'T41'],
        ['T96', 'T97', 'T98', 'T101', 'T119', 'T101'],
        ['T125', 'T46', 'T92', 'T93', 'T83'],
        ['T125', 'T54', 'T93', 'T83'],
        ['T96', 'T97', 'T98', 'T109', 'T111', 'T112', 'T113', 'T120', 'T100', 'T100'],
        ['T125', 'T67', 'T85'],
        ['T125', 'T46', 'T49', 'T50', 'T81'],
        ['T125', 'T3', 'T87', 'T88', 'T24', 'T25', 'T26', 'T22', 'T23', 'T78'],
        ['T125', 'T1', 'T2', 'T75'],
        ['T96', 'T97', 'T98', 'T99', 'T101', 'T119', 'T117', 'T118'],
        ['T125', 'T56', 'T57', 'T58', 'T83'],
        ['T125', 'T19', 'T22', 'T23', 'T78'],
        ['T96', 'T97', 'T98', 'T109', 'T120', 'T109'],
        ['T96', 'T97', 'T98', 'T99', 'T115', 'T116', 'T105', 'T123', 'T124', 'T119', 'T99'],
        ['T125', 'T67', 'T68', 'T69', 'T72', 'T73', 'T74', 'T85'],
        ['T96', 'T97', 'T98', 'T117', 'T118'],
        ['T96', 'T97', 'T98', 'T100', 'T99', 'T99'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T105', 'T106', 'T107', 'T119', 'T117', 'T118'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T102', 'T103', 'T119', 'T109', 'T111', 'T121', 'T122'],
        ['T96', 'T97', 'T98', 'T99', 'T100', 'T99'],
        ['T96', 'T97', 'T98', 'T100', 'T109', 'T111', 'T112', 'T113', 'T120', 'T100'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T119', 'T117', 'T118'],
        ['T125', 'T1', 'T86', 'T6', 'T7', 'T6'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T105', 'T123', 'T124', 'T102', 'T104', 'T119', 'T117'],
        ['T125', 'T3', 'T87', 'T16', 'T17', 'T18', 'T14', 'T15', 'T12', 'T13', 'T77'],
        ['T96', 'T97', 'T98', 'T109', 'T111', 'T112', 'T114', 'T120', 'T109'],
        ['T125', 'T27', 'T28', 'T29', 'T79'],
        ['T125', 'T38', 'T39', 'T40', 'T91', 'T47', 'T48', 'T51', 'T52', 'T53', 'T92', 'T93'],
        ['T125', 'T54', 'T82']]
        print ("pop is： ")
        for sub in pop:
            print (sub)
        generate_seq_ga_check.Elitist_clear()
        recordFun.recordInitalCaseTime()
        recordFun.recordInitalElitistInfo()
        testdata, M, V = execute.gdata(pop, SM)  # 为序列产生数据，执行测试用例，获取评估矩阵M
        SVTCM = generate_seq_ga_check.SVTCM_ICM(pop, M, 1, [])  # webtestror‘s part
        ######################## 以上为评估初始种群############################
        iteration = 0  # 迭代次数
        cflag = 0
        while True:
            print "iteration time:", iteration
            flag, cover_flag = execute.table_handle(M)
            #print ("1敏感路径入口点覆盖标识cover_flag = ", cover_flag)
            #print ("1敏感路径完全覆盖标识flag = ", flag)
            coverage, Uncover_path_list = execute.table_coverage(M)  # M是评估矩阵，基于M获取违背覆盖的敏感路径集Uncover——path——list）
            pop_time1 = datetime.now()
            pop_time = pop_time1 - pop_time0
            recordFun.recordPopTime(iteration, pop_time)  # 更新种群执行时间
            recordFun.recordPathCoverrate(iteration, coverage)  # 更新种群覆盖率
            recordFun.recordPathCover_un(iteration, Uncover_path_list)  # 更新种群未覆盖路径
            #recordFun.recordAllM(iteration, M)  # 更新评估矩阵
            #recordFun.recordPopTestdata(iteration, testdata, V)  # 更新数据和变量
            pop_time0 = datetime.now()
            ###原始####
            # if iteration > Max or cflag == 1 :
            #    if iteration > Max:
            #        recordFun.recordInitalPop(iteration, pop)  # 更新种群
            #        print "未生成完整测试套件"
            #    break
            recordFun.recordInitalPop(iteration, pop)  # 更新种群
            ###耦合####
            if iteration > Max or flag == 1:
            #if iteration > Max or cflag == 1:
                generate_seq_ga_check.recordElitist(iteration)
                recordFun.recordAllM(iteration,M)  # 更新评估矩阵
                recordFun.recordAllM(iteration,pop)
                break
            #if cover_flag == 1:  # 如果覆盖了所有敏感路径的入口点
            #    if flag == 1:  # 在判断是否覆盖所有的敏感路径
            #        testcase = zip(pop, testdata)  # 测试用例等于序列+数据
            #        break
            #    else:  # 进行局部搜索
            #        print ("第", iteration, "代" + "进行局部搜索")
            #        print ("进入局部搜索的pop", pop)
            #        print ("进入局部搜索的testdata", testdata)
            #        print ("进入局部搜索的M", M)
            #        print ("进入局部搜索的Uncover_path_list", Uncover_path_list)
            #        overgame, testdata = local_search_alg.Local_search_for_HC(pop, testdata, M, Uncover_path_list, SM)
            #        if overgame == 1:
            #            testcase = zip(pop, testdata)  # 测试用例等于序列+数据
            #            print (testcase)
            #            break  # 退出while ，执行下面的时间记录代码
            #        else:
            #            tip = "搜索失败，结束"
            #            print (tip)
            #            return tip  # 退出整个程序，不会执行下面的时间记录代码了
            #else:
            print ("继续全局GA，搜序列")
            pop,V,M,testdata,SVTCM = generate_seq_ga_check.GA_Couple(iteration, pop, pc, pm, popsize, SM, M,testdata,V,SVTCM)#webtestor
            #pop, V, M, testdata, cflag = generate_seq_ga_check.GA(iteration, pop, pc, pm, popsize, SM, M, testdata,V)  # 原始串行GA
            iteration += 1
        endTime = datetime.now()
        usertime = endTime - startTime
        print usertime
        recordFun.recordPopTime(iteration, usertime)
        LunHui += 1
        seq_to_script.clear_driver()


if __name__ == '__main__':  # not execute when import as a module
    main()

    # table = []
    # print (len(table))
    # t = [ 2, 0.0, 1.004950495049505, 0.0, 0.8544139861204133]
    # for i in range(2):
    #     table.append(t)
    # print (table)
"""
pop=[['T125', 'T38', 'T91', 'T47', 'T48', 'T51', 'T52', 'T53', 'T49', 'T50', 'T92', 'T93'],
        ['T125', 'T46', 'T92', 'T82'],
        ['T96', 'T97', 'T98', 'T109', 'T111', 'T112', 'T114', 'T110', 'T120', 'T109'],
        ['T96', 'T97', 'T98', 'T109', 'T110', 'T110'],
        ['T125', 'T38', 'T91', 'T81'],
        ['T125', 'T19', 'T89', 'T32', 'T33', 'T34', 'T35', 'T37', 'T28', 'T29', 'T30', 'T31'],
        ['T96', 'T97', 'T98', 'T109', 'T111', 'T121', 'T122', 'T110', 'T120', 'T117', 'T118'],
        ['T96', 'T97', 'T98', 'T100', 'T101', 'T102', 'T103', 'T102'],
        ['T125', 'T3', 'T8', 'T9', 'T10', 'T76'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T102', 'T104', 'T102'],
        ['T125', 'T11', 'T16', 'T17', 'T18', 'T77'],
        ['T125', 'T46', 'T81'],
        ['T125', 'T54', 'T93', 'T57', 'T58', 'T94', 'T60', 'T61', 'T95', 'T85'],
        ['T125', 'T11', 'T77'],
        ['T125', 'T56', 'T57', 'T58', 'T57'],
        ['T96', 'T97', 'T98', 'T99', 'T109', 'T111', 'T112', 'T114', 'T120', 'T101', 'T102', 'T104'],
        ['T96', 'T97', 'T98', 'T117', 'T118'],
        ['T125', 'T54', 'T93', 'T94', 'T95', 'T68', 'T69', 'T72', 'T73', 'T74', 'T70', 'T71'],
        ['T125', 'T27', 'T28', 'T29', 'T32', 'T33', 'T37', 'T32'],
        ['T125', 'T3', 'T6', 'T7', 'T8', 'T9', 'T10', 'T87', 'T77'],
        ['T125', 'T38', 'T41', 'T42', 'T91', 'T49', 'T50', 'T47', 'T48', 'T92', 'T82'],
        ['T96', 'T97', 'T98', 'T99', 'T101', 'T105', 'T106', 'T107', 'T102', 'T103', 'T105'],
        ['T96', 'T97', 'T98', 'T99', 'T115', 'T116', 'T102', 'T103', 'T105', 'T123', 'T124', 'T105'],
        ['T125', 'T38', 'T41', 'T42', 'T41'],
        ['T96', 'T97', 'T98', 'T115', 'T116', 'T119', 'T117', 'T118'],
        ['T125', 'T56', 'T94', 'T95', 'T85'],
        ['T125', 'T38', 'T41', 'T42', 'T80'],
        ['T96', 'T97', 'T98', 'T101', 'T105', 'T106', 'T107', 'T119', 'T99', 'T109', 'T111', 'T112'],
        ['T125', 'T67', 'T68', 'T69', 'T85'],
        ['T125', 'T3', 'T87', 'T12', 'T13', 'T14', 'T15', 'T77']]
"""