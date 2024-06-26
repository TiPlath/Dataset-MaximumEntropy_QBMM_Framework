# guessed parameters to be used as initial values
a1 = 10
b= -1
c= -1

ffit2(t,h)=a1*(((t-7.5)/7.5)**(b))*((h/0.00002925)**c)

# 2nd attempt
fit ffit2(x, y) 'aggregationRateVsD32VsMRT.txt' using 1:2:3 via a1, b, c

# Plot the data and the fitted function
splot 'aggregationRateVsD32VsMRT.txt' using 1:2:3 with points, ffit2(x,y)
pause -1
