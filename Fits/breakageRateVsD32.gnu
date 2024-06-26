# guessed parameters to be used as initial values
t0 = 1
b = -1
a1 = 10

ffit(d32Prime)=a1*(d32Prime)**b

fit ffit(x) 'breakageRateVsD32.txt' via a1,b

plot 'breakageRateVsD32.txt' w p, ffit(x)
pause -1