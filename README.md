# DTester
> DTester is diversity-driven web testor.  It uses the three-dimensional weight graph to organically combine various diversity, and utilizes the graph to adaptively guide the diverse generation of web test cases.

## Experimental Environment


## web application operating environment

>Mysql5.5.62

>apache2.4.39

>Please ensure the normal use of the website (Application Under Test)

## program operating environment

>python 2.7(python3 is not supported)

>Mozilla Firefox 64.0

>geckodriver-v0.24.0


## Execution


`handle` contains the main executable files

`module` stores the Client Behavior Model of the application under test

`dataset` stores the data that needs to be collected when testing the application under test. You can modify $recordFun.py$ under handle according to your needs.

`spath` stores vulnerable path file of application under test for verification


Please ensure that the path and config to the application under test is correctly set in the code. Roughly need to modify the following files：`config.py`，`case.py`，`seq_to_script.py`

******
## phpaaCMS
data [link](https://www.aliyundrive.com/s/hDnut5RTSXm "phpaaCMS") 
