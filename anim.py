import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

dx = 0.1
datalocation = "../Case 1/Plane Output/"

m = 4.0 + 1e-10
n = 3.0 + 1e-10

for i in range(1,201):

    j = i*10
    if i==0: j=1
    # j = i
    print('%d / %d' %(i,600))

    filename = datalocation + 'planesoutput-' + str(j).zfill(4)
    table = pd.read_csv(filename, delimiter=',')
    array = np.array(table)

    x = array[:,1]
    z = array[:,3]
    xvelocity = array[:,8]
    yvelocity = array[:,9]
    zvelocity = array[:,10]
    temperature = array[:,11]
    dpm_concentration = array[:,13]
    prob_infection = array[:,20]

    grid_x, grid_z = np.mgrid[dx/2:m-dx/2:dx, dx/2:n-dx/2:dx]

    grid_u = griddata((x, z), xvelocity, (grid_x, grid_z), method='cubic')
    grid_v = griddata((x, z), yvelocity, (grid_x, grid_z), method='cubic')
    grid_w = griddata((x, z), zvelocity, (grid_x, grid_z), method='cubic')
    grid_t = griddata((x, z), temperature, (grid_x, grid_z), method='cubic')
    grid_d = griddata((x, z), dpm_concentration, (grid_x, grid_z), method='cubic')
    grid_p = griddata((x, z), prob_infection, (grid_x, grid_z), method='cubic')

    grid_u[np.isnan(grid_u)] = 0
    grid_v[np.isnan(grid_v)] = 0
    grid_w[np.isnan(grid_w)] = 0
    grid_t[np.isnan(grid_t)] = 294
    grid_d[np.isnan(grid_d)] = 0
    grid_p[np.isnan(grid_p)] = 0

    grid_d = grid_d + 1e-10
    grid_d = np.log10(np.abs(grid_d))
    grid_p.clip(min=0)

    fig,axs = plt.subplots(2,3,num = 1, figsize = (18,9), dpi=150)
    axs = axs.flatten() 

    im1 = axs[0].imshow(grid_u.T,extent=[0,m,0,n],vmin=-0.2, vmax=0.2, cmap='RdBu_r', origin='lower')
    axs[0].set_xlabel('X [m]')
    axs[0].set_ylabel('Z [m]')
    axs[0].set_title('X - Velocity [m/s]')
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im1, cax=cax)

    im2 = axs[1].imshow(grid_v.T,extent=[0,m,0,n],vmin=-0.3, vmax=0.3, cmap='RdBu_r',origin='lower')
    axs[1].set_xlabel('X [m]')
    axs[1].set_ylabel('Z [m]')
    axs[1].set_title('Y - Velocity [m/s]')
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im2, cax=cax)

    im3 = axs[2].imshow(grid_w.T,extent=[0,m,0,n],vmin=-0.6, vmax=0.6, cmap='RdBu_r',origin='lower')
    axs[2].set_xlabel('X [m]')
    axs[2].set_ylabel('Z [m]')
    axs[2].set_title('Z - Velocity [m/s]')
    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im3, cax=cax)

    im4 = axs[3].imshow(grid_t.T,extent=[0,m,0,n], vmin=294, vmax=296, cmap='Oranges',origin='lower')
    axs[3].set_xlabel('X [m]')
    axs[3].set_ylabel('Z [m]')
    axs[3].set_title('Temperature [K]')
    divider = make_axes_locatable(axs[3])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im4, cax=cax)

    im5 = axs[4].imshow(grid_d.T,extent=[0,m,0,n], vmin=-9, vmax=-5.5, cmap='Oranges',origin='lower')
    axs[4].set_xlabel('X [m]')
    axs[4].set_ylabel('Z [m]')
    axs[4].set_title('log10(DPM Concentration) [log10(kg/m^3)]')
    divider = make_axes_locatable(axs[4])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im5, cax=cax)

    im6 = axs[5].imshow(grid_p.T,extent=[0,m,0,n], vmin=0, vmax=5, cmap='Oranges',origin='lower')
    axs[5].set_xlabel('X [m]')
    axs[5].set_ylabel('Z [m]')
    axs[5].set_title('Probability of Infection [%]')
    divider = make_axes_locatable(axs[5])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im6, cax=cax)

    plt.suptitle('t=%.2f s' %(j/10), fontsize = 20)
    plt.tight_layout()
    plt.savefig('../Case 0/Plots/%d_log.png' %j)
    
    plt.close(fig)


    
    

