import numpy as np
import matplotlib.pyplot as plt

data = []
sub_data = []

with open("outfile.txt", 'r') as file:
    for line in file:
        if line[0] == "-":
            data.append(sub_data)
            sub_data.clear
        else:
            sub_data.append(line.split(","))

#print(data[0])
xmax = 0
ymax = 0
for coord in data[0]:
    if float(coord[0]) > xmax:
        xmax=float(coord[0])
    if float(coord[1]) > ymax:
        ymax=float(coord[1])

print(xmax)
print(ymax)

#next step is to correlate x's and y's to indices
            
data = np.random.rand(10,10)

plt.imshow(data, cmap='binary')

plt.colorbar()

plt.show()