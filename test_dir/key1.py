import sys,re,os

countggg=2161
with open(sys.argv[1]) as key :
    for i in key :
        if countggg in int(i.split(",")[1]) :
            print(i)
