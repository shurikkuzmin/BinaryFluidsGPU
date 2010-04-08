set term gif
set output "$0.gif"
set zrange [2.0:5.0]
splot [0:255] [0:255] '$0' matrix w l
#!kuickshow test.jpg
reset