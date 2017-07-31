import matplotlib.pyplot as plt
import numpy as np


num_solutions = 5

energy = []
solutions = []
for i in range(num_solutions):
  solutions.append([])

input_dir = ''
file_name = 'results.txt'
input_file = open(input_dir+file_name, 'r')

for line in input_file:
    data = line.split()
    for i in range(len(data)):
        data[i] = float(data[i])

    energy.append(data[0])
    for i in range(0, num_solutions):
        solutions[i].append(data[i+1])


for solution in solutions:
    plt.loglog(energy, solution)
    plt.xlabel('Energy [normalised]')
    plt.ylabel('Solution [normalised]')
plt.axis([0.01, 30, 1e-3, 1e3])
plt.savefig('t_LOD')
plt.show()
