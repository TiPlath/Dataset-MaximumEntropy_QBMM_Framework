# guessed parameters to be used as initial values
a1 = 10
b= -1
# a1 is the rate, t0 time in which the particle number halves
ffit2(t)=a1*((t-7.5)/7.5)**(b)
#plot 'aggregationRateVsMRT.txt' w p, ffit(x)
#pause -1
#
# 2nd attempt
fit ffit2(x) 'aggregationRateVsMRT.txt' via a1,b
plot 'aggregationRateVsMRT.txt' w p, ffit2(x)
pause -1
