# guessed parameters to be used as initial values
a0=0.1 # zero is not allowed!
b= 1
c= 1
#
ffit(t)=a0 + t/b
plot 'd32vsLS.txt' w p, ffit(x)
pause -1
#
# 2nd attempt
fit ffit(x) 'd32vsLS.txt' via a0,b
plot 'd32vsLS.txt' w p, ffit(x)
pause -1