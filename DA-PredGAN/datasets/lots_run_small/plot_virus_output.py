import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt

import numpy as np
import sys, os


filename = sys.argv[1]
filebase = os.path.splitext(filename)[0]
print filebase
#k = sys.argv[1]
#filename = 'group-output-time' + str(k) + '.csv'


f = open(filename, "r")
contents = f.readlines()
f.close()

count_row = 0
ngroups = 0 
for iline,line in enumerate(contents):
    line_contents = line.split()
    first_word = line_contents[0]
    if first_word[0:1] == '#':
        continue
    elif first_word[0:3] == 'GRO':
        ngroups = ngroups + 1
        count_row = 0 
        continue
    elif first_word[0:3] == 'min' or first_word[0:3] == 'sum':
        continue

    count_col = len(line_contents)
    #people_density[count,:,ngroups-1] = np.array(line_contents,dtype=np.float64)
    count_row = count_row + 1

print ('nx by ny:', count_row, count_col)
print('ngroups:',ngroups)

people_density = np.zeros((count_row, count_col, ngroups))

ngroups = 0
count_row = 0
for iline,line in enumerate(contents):
    line_contents = line.split()
    first_word = line_contents[0]
    if first_word[0:1] == '#':
        if line_contents[1][0:4] == 'TIME':
            acc_time = np.float(line_contents[2])
        continue
    elif first_word[0:3] == 'GRO':
        ngroups = ngroups + 1
        count_row = 0
        continue
    elif first_word[0:3] == 'min' or first_word[0:3] == 'sum':
        continue 
    people_density[count_row,:,ngroups-1] = np.array(line_contents,dtype=np.float64)
    count_row = count_row + 1


#plt.figure(figsize=(20,10))

#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)

fig, axs = plt.subplots(2, 4, figsize=(20,10))
contour_levels = 10

seconds_in_day = 3600*24.
time = "{:.2f}".format(acc_time / seconds_in_day)

fig.suptitle('Day ' + str(time) )
#cs = ax.contourf(X, Y, Z, levels, cmap=cmap, extend=extend, origin=origin)
#fig.colorbar(cs, ax=ax, shrink=0.9)

cs = axs[0, 0].contourf(people_density[:,:,0],contour_levels) # plot(x, y)
axs[0, 0].set_title('Home - S')
axs[0, 0].set_aspect('equal')
fig.colorbar(cs, ax=axs[0, 0], shrink=0.7) # orientation='horizontal', 
print('0,0, min and max', np.min(people_density[:,:,0]), np.max(people_density[:,:,0]))

cs = axs[0, 1].contourf(people_density[:,:,1],contour_levels)#plot(x, y, 'tab:orange')
axs[0, 1].set_title('Home - E')
axs[0, 1].set_aspect('equal')
fig.colorbar(cs, ax=axs[0, 1], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,1]), np.max(people_density[:,:,1]))

cs = axs[0, 2].contourf(people_density[:,:,2],contour_levels)#plot(x, y, 'tab:orange')
axs[0, 2].set_title('Home - I')
axs[0, 2].set_aspect('equal')
fig.colorbar(cs, ax=axs[0, 2], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,2]), np.max(people_density[:,:,2]))

cs = axs[0, 3].contourf(people_density[:,:,3],contour_levels)#plot(x, y, 'tab:orange')
axs[0, 3].set_title('Home - R')
axs[0, 3].set_aspect('equal')
fig.colorbar(cs, ax=axs[0, 3], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,3]), np.max(people_density[:,:,3]))

cs = axs[1, 0].contourf(people_density[:,:,4],contour_levels)#plot(x, -y, 'tab:green')
axs[1, 0].set_title('Mobile - S')
axs[1, 0].set_aspect('equal')
fig.colorbar(cs, ax=axs[1, 0], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,4]), np.max(people_density[:,:,4]))

cs = axs[1, 1].contourf(people_density[:,:,5],contour_levels)#plot(x, -y, 'tab:red')
axs[1, 1].set_title('Mobile - E')
axs[1, 1].set_aspect('equal')
fig.colorbar(cs, ax=axs[1, 1], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,5]), np.max(people_density[:,:,5]))

cs = axs[1, 2].contourf(people_density[:,:,6],contour_levels)#plot(x, -y, 'tab:red')
axs[1, 2].set_title('Mobile - I')
axs[1, 2].set_aspect('equal')
fig.colorbar(cs, ax=axs[1, 2], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,6]), np.max(people_density[:,:,6]))

cs = axs[1, 3].contourf(people_density[:,:,7],contour_levels)#plot(x, y, 'tab:orange')
axs[1, 3].set_title('Mobile - R')
axs[1, 3].set_aspect('equal')
fig.colorbar(cs, ax=axs[1, 3], shrink=0.7)
print('0,0, min and max', np.min(people_density[:,:,7]), np.max(people_density[:,:,7]))

#for ax in axs.flat:
#    ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
#plt.contourf(people_density_gp1,10)
#plt.show()
#plt.savefig('virus-figure-day-' + k + '.png')filename
plt.savefig(filebase + '.png')






