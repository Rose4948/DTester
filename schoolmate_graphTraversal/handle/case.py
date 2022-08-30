# coding=utf-8
'''定义各种测试用例'''
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
import random
import time
from selenium.webdriver.common.keys import Keys
#schoolmate的登录
def dologin(driver, user, psw): # 登录，参数是用户名和密码  T3
    driver.find_element_by_name("username").clear()
    driver.find_element_by_name("username").send_keys(user)
    # element is hidden ,must execute js to obtain
    # driver.execute_script("document.getElementById('input-login-password').style.display='block';")
    driver.find_element_by_name("password").clear()
    driver.find_element_by_name("password").send_keys(psw)
    driver.find_element_by_xpath("html/body/table[2]/tbody/tr[2]/td[3]/table/tbody/tr/td/table/tbody/tr[1]/td[2]/div/form/table/tbody/tr[3]/td[2]/input").click()


def doquit(driver):  # 退出火狐
    driver.close()
    driver.quit()

def T19(driver,object):
    driver.find_element_by_xpath(object).click()
    Select(driver.find_element_by_name("semester")).select_by_visible_text("All")
def T24(driver,object):
    # s = driver.find_element_by_name("semester")
    # s.click()
    # Select(s).select_by_visible_text("All")
    Select(driver.find_element_by_name("semester")).select_by_visible_text("All")
    driver.find_element_by_name("delete[]").click()
    time.sleep(3)
    # driver.find_element_by_xpath(object).click()
###################################################################3333
def click_event_id(driver, click_object):
    driver.find_element_by_id(click_object).click()

def click_event_xpath(driver, click_object):
    driver.find_element_by_xpath(click_object).click()
    
def click_event_name(driver, click_object):
    driver.find_element_by_name(click_object).click()

def click_event_link_all(driver,click_object):
    driver.find_element_by_link_text(click_object).click()


def click_event_css(driver, click_object):
    driver.find_element_by_css_selector(click_object).click()

################  以上为基础封装  ###########################################

def click_event_input_one(driver,p1,v1,click_object):
    # driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_xpath(click_object).click()


def click_event_input_two(driver,p1,p2,v1,v2,click_object):
    # driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    # driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_xpath(click_object).click()


def click_event_input_three(driver,p1,p2,p3,v1,v2,v3,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_xpath(click_object).click()

def click_event_input_four_add(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).send_keys(v3)

    # driver.find_element_by_name(p4).send_keys(v4)
    # randlist = [1,3,4]
    # rand = random.choice(randlist)
    Select(driver.find_element_by_name(p4)).select_by_index(3)

    driver.find_element_by_xpath(click_object).click()
def click_event_input_four_parent(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)

    Select(driver.find_element_by_name(p3)).select_by_index(0)
    Select(driver.find_element_by_name(p4)).select_by_index(0)

    driver.find_element_by_xpath(click_object).click()


def click_event_input_four(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).send_keys(v3)

    # driver.find_element_by_name(p4).send_keys(v4)
    Select(driver.find_element_by_name(p4)).select_by_index(0)

    driver.find_element_by_xpath(click_object).click()


def click_event_input_four_teacher(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_xpath(click_object).click()

def click_event_input_four_faq(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    Select(driver.find_element_by_name(p3)).select_by_index(0)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_xpath(click_object).click()


def click_event_input_four_webchess(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).send_keys(v2)
    Select(driver.find_element_by_name(p3)).select_by_index(0)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)

    driver.find_element_by_xpath(click_object).click()

def click_event_input_five(driver,p1,p2,p3,p4,p5,v1,v2,v3,v4,v5,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    driver.find_element_by_xpath(click_object).click()


def click_event_input_five_faq(driver,p1,p2,p3,p4,p5,v1,v2,v3,v4,v5,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    # driver.find_element_by_name(p3).send_keys(v3)
    Select(driver.find_element_by_name(p3)).select_by_index(0)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).send_keys(v5)
    driver.find_element_by_xpath(click_object).click()

def click_event_input_six(driver,p1,p2,p3,p4,p5,p6,v1,v2,v3,v4,v5,v6,click_object):
    # driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    # driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    # driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    # driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    # driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    # driver.find_element_by_name(p6).clear()
    driver.find_element_by_name(p6).send_keys(v6)
    driver.find_element_by_xpath(click_object).click()


def click_event_input_sevent(driver,p1,p2,p3,p4,p5,p6,p7,v1,v2,v3,v4,v5,v6,v7,click_object):
    # driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    # driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    # driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    # driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    # driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    # driver.find_element_by_name(p6).clear()
    driver.find_element_by_name(p6).send_keys(v6)
    driver.find_element_by_name(p7).send_keys(v7)
    driver.find_element_by_xpath(click_object).click()


def webchess_T3(driver,p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v4)

    driver.find_element_by_name(p6).click()  #choose checkbox

    driver.find_element_by_name(p7).click()

    driver.find_element_by_name(p8).click()

    driver.find_element_by_name(p9).clear()
    driver.find_element_by_name(p9).send_keys(v9)
    driver.find_element_by_name(click_object).click()


def webchess_T4(driver,p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)

    driver.find_element_by_name(p6).click()  #choose checkbox

    driver.find_element_by_name(p7).click()

    driver.find_element_by_name(p8).click()

    driver.find_element_by_name(p9).clear()
    driver.find_element_by_name(p9).send_keys(v9)
    driver.find_element_by_name(click_object).click()
    alert = driver.switch_to_alert()
    alert.accept()



def click_event_input_12(driver,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    driver.find_element_by_name(p6).clear()
    driver.find_element_by_name(p6).send_keys(v6)
    driver.find_element_by_name(p7).clear()
    driver.find_element_by_name(p7).send_keys(v7)
    driver.find_element_by_name(p8).clear()
    driver.find_element_by_name(p8).send_keys(v8)
    driver.find_element_by_name(p9).clear()
    driver.find_element_by_name(p9).send_keys(v9)
    driver.find_element_by_name(p10).clear()
    driver.find_element_by_name(p10).send_keys(v10)
    driver.find_element_by_name(p11).clear()
    driver.find_element_by_name(p11).send_keys(v11)
    driver.find_element_by_name(p12).clear()
    driver.find_element_by_name(p12).send_keys(v12)
    driver.find_element_by_xpath(click_object).click()

#论文示例
def onclick_event_input_two(driver,C1,C2,data1,data2,C):
    # driver.find_element_by_name(C1).clear()
    driver.find_element_by_name(C1).send_keys(data1)
    # driver.find_element_by_name(C2).clear()
    driver.find_element_by_name(C2).send_keys(data2)
    driver.find_element_by_xpath(C).click()
#schoolmate的点击弹框link
def link_text(driver):
    driver.find_element(By.LINK_TEXT, "LCbr3g").click()

def click_link(driver):
    # driver.find_element_by_link_text(click_object).click()
    alert = driver.switch_to_alert()
    alert.accept()
#schoolmate的点击弹框link
def click_event_links(driver):
    alert = driver.switch_to_alert()
    time.sleep(2)
    alert.accept()

#下拉输入
def down_click_three(driver,p1,p2,p3,click_object):
    olist = []
    student = driver.find_element_by_name("student")
    student.click() #点击student选择框
    Select(student).select_by_index(2)
    Select(driver.find_element_by_name(p2)).select_by_index(0)
    semeter = driver.find_element_by_name(p2).find_element_by_tag_name("option")
    while 1:
        select = driver.find_element_by_name(p3)
        olist = select.find_element_by_tag_name("option")
        if olist :
            Select(driver.find_element_by_name(p3)).select_by_index(0)
            break
        else:
            r = random.randint(semeter)
            Select(driver.find_element_by_name(p2)).select_by_index(r)
        break
    driver.find_element_by_xpath(click_object).click()


def down_click_two(driver,p1,p2,click_object):
    Select(driver.find_element_by_name(p1)).select_by_index(0)
    Select(driver.find_element_by_name(p2)).select_by_index(0)
    driver.find_element_by_xpath(click_object).click()


def down_click_xpath_one(driver,p1,click_object):
    Select(driver.find_element_by_name(p1)).select_by_index(1)
    driver.find_element_by_xpath(click_object).click()


def down_click_name_one(driver,p1,click_object):
    Select(driver.find_element_by_name(p1)).select_by_index(0)
    driver.find_element_by_name(click_object).click()

def down_click_name_one_T97(driver,p1,click_object):
    print "#T97"
    driver.find_element_by_name(p1).click()
    dropdown = driver.find_element(By.NAME, "semester")
    dropdown.find_element(By.XPATH, "//option[. = 's5UBs5']").click()

def click_event_input_css_five(driver,p1,p2,p3,p4,p5,v1,v2,v3,v4,v5,click_object):
    driver.find_element_by_name(p1).send_keys(v1)
    # driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    # driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    # driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    # driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    driver.find_element_by_css_selector(click_object).click()

def click_event_input_name_four(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).send_keys(v1)
    # driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    # driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    # driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    # driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(click_object).click()

def faqfore_T2_T3(driver,p1,p2,v1,v2,click_object):
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(click_object).click()


def webchess_T1(driver,p1,p2,v1,v2,click_object):
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(click_object).click()


def webchess_T19(driver,p1,p2,v1,v2,click_object):
    # Select(driver.find_element_by_name(p1)).select_by_visible_text("test")
    Select(driver.find_element_by_name(p1)).select_by_index(2)
    time.sleep(1)
    driver.find_element_by_xpath("(//input[@name='color'])[2]").click()
    driver.find_element_by_xpath(click_object).click()

def webchess_T18(driver,p1,v1,click_object):
    sreach_windows = driver.current_window_handle
    Select(driver.find_element_by_name(p1)).select_by_index(1)
    time.sleep(1)
    driver.find_element_by_xpath(click_object).click()
    time.sleep(3)
    all_handles = driver.window_handles
    return sreach_windows,all_handles

def webchess_T21(driver,p1,p2,v1,v2,click_object,window,all_window):
    for handle in all_window:
        if handle != window:
            driver.switch_to_window(handle)
            driver.find_element_by_name(p1).send_keys(v1)
            driver.find_element_by_name(p2).send_keys(v2)
            driver.find_element_by_name(click_object).click()
            break
    driver.switch_to_window(all_window[0])
    time.sleep(1)
def webchess_T17(driver,p1,p2,p3,p4,p5,v1,v2,v3,v4,v5,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(1)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(1)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(1)
    driver.find_element_by_xpath(click_object).click()

#addressbook的登录
def login(driver, user, psw): # 登录，参数是用户名和密码  T3
    driver.find_element_by_name("user").clear()
    driver.find_element_by_name("user").send_keys(user)
    # element is hidden ,must execute js to obtain
    #driver.execute_script("document.getElementById('input-login-password').style.display='block';")
    driver.find_element_by_name("pass").clear()
    driver.find_element_by_name("pass").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div/div[4]/form/input[3]").click()
#openconf的登录
def loginOpen(driver): # 登录，参数是用户名和密码  T3
    driver.find_element_by_name("uname").clear()
    driver.find_element_by_name("uname").send_keys(user)
    # element is hidden ,must execute js to obtain
    #driver.execute_script("document.getElementById('input-login-password').style.display='block';")
    driver.find_element_by_name("upwd").clear()
    driver.find_element_by_name("upwd").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div/div[4]/form/input[3]").click()

def addressbook_Ninput(driver,t,induval,j,click_object):
    sp = ["bday","bmonth","aday","amonth","new_group"]
    for i in range(len(t)):
        if t[i] == "photo":   #上传文件
            driver.find_element_by_name("photo").clear()
            driver.find_element_by_name("photo").send_keys("/home/wushumei/python/graphTraversal_addressbook/support/upload.csv")
        elif t[i] in sp:
            Select(driver.find_element_by_name(t[i])).select_by_index(4)
        else:
            driver.find_element_by_name(t[i]).clear()
            driver.find_element_by_name(t[i]).send_keys(induval[j-len(t)+i])

    driver.find_element_by_xpath(click_object).click()

def upaload(driver,p1,click_object):#上传操作
    driver.find_element_by_id(p1).clear()
    driver.find_element_by_id(p1).send_keys("/home/wushumei/python/graphTraversal_addressbook/support/upload.csv")
    driver.find_element_by_name(click_object).click()


def addressbook_four_input(driver,p1,p2,p3,p4,v1,v2,v3,v4,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    Select(driver.find_element_by_name(p2)).select_by_index(1)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v2)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(click_object).click()


def click_name_one_input(driver,p1,v1,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_xpath(click_object).click()

#phpcss 特殊处理case
def loginForPhpCss(driver,username,password):
    driver.find_element_by_name("username").clear()
    driver.find_element_by_name("username").send_keys(username)
    driver.find_element_by_name("password").clear()
    driver.find_element_by_name("password").send_keys(password)
    driver.find_element_by_name("image").click()


def php_upload(driver,p1,url,click_object):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(url)
    driver.find_element_by_name(click_object).click()


def php_name_one(driver,p1,v1,click_obejct):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(click_obejct).click()


def php_name_two(driver,p1,p2,v1,v2,click_obejct):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(click_obejct).click()


def php_name_three(driver,p1,p2,p3,v1,v2,v3,click_obejct):
    driver.find_element_by_name(p1).send_keys(v1)

    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(click_obejct).click()


def php_name_six(driver,p1,p2,p3,p4,p5,p6,v1,v2,v3,v4,v5,v6,click_obejct):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    driver.find_element_by_name(p3).clear()
    driver.find_element_by_name(p3).send_keys(v3)
    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    driver.find_element_by_name(p6).clear()
    driver.find_element_by_name(p6).send_keys(v6)
    driver.find_element_by_name(click_obejct).click()
    alert = driver.switch_to_alert()
    alert.accept()
    # driver.find_element_by_link_text("OK").click()


def php_name_eight(driver,p1,p2,p3,p4,p5,p6,p7,p8,v1,v2,v3,v4,v5,v6,v7,v8,click_obejct):
    driver.find_element_by_name(p1).clear()
    driver.find_element_by_name(p1).send_keys(v1)
    driver.find_element_by_name(p2).clear()
    driver.find_element_by_name(p2).send_keys(v2)
    #上传文件
    driver.find_element_by_id(p3).clear()
    driver.find_element_by_id(p3).send_keys("/home/wushumei/python/graphTraversal_addressbook/support/888.png")

    driver.find_element_by_name(p4).clear()
    driver.find_element_by_name(p4).send_keys(v4)
    driver.find_element_by_name(p5).clear()
    driver.find_element_by_name(p5).send_keys(v5)
    #选择下拉框
    s=driver.find_element_by_name(p6) #定位到下拉框
    s.find_element_by_xpath("/html/body/form/table[2]/tbody/tr[5]/td[2]/select/option[2]").click()

    driver.find_element_by_name(p7).clear()
    driver.find_element_by_name(p7).send_keys(v7)
    driver.find_element_by_name(p8).clear()
    driver.find_element_by_name(p8).send_keys(v8)
    driver.find_element_by_name(click_obejct).click()


def php_id_one(driver,p1,v1,click_object):
    driver.find_element_by_id(p1).clear()
    driver.find_element_by_id(p1).send_keys(v1)
    driver.find_element_by_xpath(click_object).click()


def php_link(driver,click_object):

    driver.find_element_by_link_text(click_object).click()



def T62(driver):
    driver.find_element_by_xpath("(//input[@name='restorefrom'])[2]").click()


def T61(driver):
    driver.find_element_by_xpath("(//input[@name='bfzl'])[2]").click()
    b = driver.find_element_by_name("tablename")
    b.find_element_by_xpath("/html/body/form/table/tbody/tr[4]/td[1]/select/option[2]").click()
    driver.find_element_by_xpath("(//input[@name='weizhi'])[2]").click()
    driver.find_element_by_name("act").click()

#*******************************************************************openconf*************************************************************************************************************************
def NO_input_click(driver, method, click_object):
    if method == "id":
        driver.find_element_by_id(click_object).click()
    elif method == "xpath":
        driver.find_element_by_xpath(click_object).click()
    elif method == "name":
        driver.find_element_by_name(click_object).click()
    elif method == "css":
        driver.find_element_by_css_selector(click_object).click()
    elif method == "link":
        driver.find_element_by_link_text(click_object).click()
#*******************两个输入的登陆处理函数**********************************************************
def openconf_Two_Login(driver, t1, t2, induval1, induval2, method, click_object):
    driver.find_element_by_name(t1).clear()
    driver.find_element_by_name(t1).send_keys(induval1)
    driver.find_element_by_name(t2).clear()
    driver.find_element_by_name(t2).send_keys(induval2)
    if method=="id":
        driver.find_element_by_id(click_object).click()
    elif method=="xpath":
        driver.find_element_by_xpath(click_object).click()
    elif method=="name":
        driver.find_element_by_name(click_object).click()
    elif method=="css":
        driver.find_element_by_css_selector(click_object).click()
    elif method=="link":
        driver.find_element_by_link_text(click_object).click()
#****************************一般输入的处理函数*************************************************
def openconf_Ninput(driver, t, induval, j, method, click_object):
    if len(t) in [6,18]:
        induval[j-1]=induval[j-2]
    elif len(t) == 9:
        induval[j - 1] = "test12345"
        induval[j - 2] = ""
        induval[j - 3] = ""
    for i in range(len(t)):
        tt = t[i].split("-")
        if len(tt) < 2:#解决了T18的问题
            if t[i]=='name_last':#输入框中有默认值，clear清除后，无法重新赋值，就用此法解决
                driver.find_element_by_id(t[i]).send_keys(Keys.CONTROL, 'a')
                driver.find_element_by_id(t[i]).send_keys(induval[j - len(t) + i])
            else:
                driver.find_element_by_name(t[i]).clear()
                driver.find_element_by_name(t[i]).send_keys(induval[j - len(t) + i])
        elif tt[1] == "click":
            driver.find_element_by_name(tt[0]).click()  #choose checkbox
        elif tt[1] == "id_click":
            if tt[0]=="OC_upload_open":
                v2=driver.find_element_by_id("OC_upload_open2").is_selected()
                if v2 == 1:
                    tt[0]="OC_upload_open1"
                else:
                    tt[0]="OC_upload_open2"
            driver.find_element_by_id(tt[0]).click()  # choose checkbox
        elif tt[1] == "select":
            if tt[0]=='recipient':
                Select(driver.find_element_by_name(tt[0])).select_by_index(1)
            else:
                Select(driver.find_element_by_name(tt[0])).select_by_index(0)
        elif tt[1] == "addSelection":
            # -addSelection//需要进行处理，10-14号再说
            Select(driver.find_element_by_name(tt[0])).select_by_index(0)
        else:
            #不需要单独考虑name-5，short-5但需要考虑paa-int
            if tt[1]=="5" :
                tt[0]=tt[0]+"-"+tt[1]
            driver.find_element_by_name(tt[0]).clear()
            driver.find_element_by_name(tt[0]).send_keys(induval[j-len(t)+i])
    if method=="id":
        driver.find_element_by_id(click_object).click()
    elif method=="xpath":
        driver.find_element_by_xpath(click_object).click()
    elif method=="name":
        driver.find_element_by_name(click_object).click()
    elif method=="css":
        driver.find_element_by_css_selector(click_object).click()
    elif method=="link":
        driver.find_element_by_link_text(click_object).click()













