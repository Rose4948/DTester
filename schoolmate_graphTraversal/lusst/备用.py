if url == "http://localhost/openconf/author/submit.php":
    pass
elif url == "http://localhost/openconf/author/edit.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("passwordfld").clear()
    driver.find_element_by_name("passwordfld").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/p/input").click()
elif url == "http://localhost/openconf/author/upload.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("pwd").clear()
    driver.find_element_by_name("pwd").send_keys(psw)
    driver.find_element_by_name("pwd").clear()
    driver.find_element_by_name("pwd").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/div/input").click()
elif url == "http://localhost/openconf/author/paper.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("passwordfld").clear()
    driver.find_element_by_name("passwordfld").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/p/input").click()
elif url == "http://localhost/openconf/author/withdraw.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("passwordfld").clear()
    driver.find_element_by_name("passwordfld").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/p/input").click()
elif url == "http://localhost/openconf/author/status.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("passwordfld").clear()
    driver.find_element_by_name("passwordfld").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/p/input").click()
elif url == "http://localhost/openconf/review/signin.php":
    driver.find_element_by_name("pid").clear()
    driver.find_element_by_name("pid").send_keys(user)
    driver.find_element_by_name("passwordfld").clear()
    driver.find_element_by_name("passwordfld").send_keys(psw)
    driver.find_element_by_xpath("/html/body/div[4]/form/p/input").click()
elif url == "http://localhost/openconf/chair/signin.php":
    driver.find_element_by_name("uname").clear()
    driver.find_element_by_name("uname").send_keys("chair")
    driver.find_element_by_name("upwd").clear()
    driver.find_element_by_name("upwd").send_keys("chair12345")
    driver.find_element_by_xpath("/html/body/div[4]/form/table/tbody/tr[3]/th/input").click()
elif url == "http://localhost/openconf/":
    driver.find_element_by_name("keycode").clear()
    driver.find_element_by_name("keycode").send_keys("revkey")
    driver.find_element_by_xpath("/html/body/div[4]/ul[2]/li[2]/form/input[2]").click()