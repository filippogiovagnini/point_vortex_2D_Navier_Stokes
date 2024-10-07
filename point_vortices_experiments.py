import matplotlib.pyplot as plt
import numpy as np

# Compute the Green's function in (x, y) and (xi, yi)
def greenfunc(x, y, xi, yi, rsmall):
    # rsmall is a small number to avoid singularity
    rsq = (x-xi)**2 + (y-yi)**2
    rsqpos = max([rsq, rsmall])
    a = 0.0005
    p = 4
    greenx = (1.0/4*np.pi)*1/p*np.log(a*rsqpos**p + 1) + (1.0/4*np.pi)*np.log(rsqpos)
    return greenx

# Compute the streamfunction at a point (x, y) due to a set of point vortices
def streamfunction(xarr, yarr, carr, rsmall, x, y):
    nv = len(carr)
    streamx = 0
    for i in range(nv):
        streamx = streamx + carr[i]*greenfunc(x, y, xarr[i], yarr[i], rsmall)
    return streamx

# Compute the velocity at a point (x, y) due to a set of point vortices
def Velocity_at_point(xarr, yarr, carr, rsmall, x, y):
    
    dy = 0.05; dx = 0.05
    u_vel = -(streamfunction(xarr, yarr, carr, rsmall, x, y+dy) - streamfunction(xarr, yarr, carr, rsmall, x, y-dy))/2*dy
    v_vel =  (streamfunction(xarr, yarr, carr, rsmall, x+dx, y) - streamfunction(xarr, yarr, carr, rsmall, x-dx, y))/2*dx
    return u_vel, v_vel

# Update the position of the point vortices using the Forward Euler method for the first time step
def FE(xarr, yarr, carr, rsmall):
    dt= 1
    for i in range(0,len(xarr)):
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        xarr[i] = xarr[i] + dt*u_vel
        yarr[i] = yarr[i] + dt*v_vel
    return xarr, yarr

# Update the position of the point vortices using the Forward Euler method for the second time step onwards
def Leapfrog(xarr, yarr,xarr_old,yarr_old, carr, rsmall):
    dt = 1
    for i in range(0,len(xarr)):
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        xarr[i] = xarr_old[i] + 2*dt*u_vel
        yarr[i] = yarr_old[i] + 2*dt*v_vel
    return xarr, yarr

# Compute the streamfield due to a set of point vortices and return a matrix, therefore not just in (x, y) as above
def Update_stream_field(xarr, yarr, carr, rsmall, xg, yg, nr):
    ng = 2*nr+1
    streamfield = np.zeros((ng, ng))
    for j in range(ng):
        for i in range(ng):
            streamx = streamfunction(xarr, yarr, carr, rsmall, xg[i], yg[j])
            streamfield[i, j] = streamx
    return streamfield


def plotfield2(xg, yg, field, fieldname, xarr, yarr):
    '''
    Plot 2-D field on (x, y) plane.
    Input: xg - grid points in x
           yg - grid points in y
           field - 2-d array of field values at grid points
           fieldname - name of field
    Output: Plot of field
    '''
    plt.title(fieldname)    
    plt.xlabel('x')
    plt.ylabel('y')

    # number of contour levels to use
    nlevs = 20

    # Note that plt.contourf plots an array as transposed due to convention on 
    # first index as columns and second as rows. So to get a
    # plot of field in expected orientation need to transpose the array (using .T).
    plt.contourf(xg, yg, field.T, nlevs)

    # This is to add the black dots for the point vortices
    plt.scatter(xarr, yarr, color='black')
    return
    
if __name__ == '__main__':
    # Create a grid for plotting the streamfunction
    nr = 10
    drscal = 0.2
    ng = 2*nr+1
    xg = drscal*(np.arange(ng)-nr)
    yg = np.copy(xg)
    streamfield = np.zeros([ng, ng])
    streamname = 'streamfunction'

    # Parameter for the Green's function
    rsmall = 1.e-10

    # Number of point vortices
    nv=11

    theta = np.linspace(0,2*np.pi,nv)
    xarr = np.cos(theta)
    yarr = np.sin(theta)

    # Strength of the point vortices
    carr = np.ones(nv) +np.random.randn(nv)

    k = 0
    nt = 10000
    xstore = np.zeros([nt,nv])
    ystore = np.zeros([nt,nv])
    xstore[0,:] = xarr
    ystore[0,:] = yarr

    plt.figure(1)

    if k==0:
        k = k+1
        xarr, yarr= FE(xarr, yarr, carr, rsmall)
        xstore[k,:] = xarr
        ystore[k,:] = yarr

    while (k < nt-1):

        xstore[k+1,:], ystore[k+1,:]= Leapfrog(xstore[k,:], ystore[k,:],xstore[k-1,:], ystore[k-1,:], carr, rsmall)
        streamfield = Update_stream_field(xstore[k+1,:],ystore[k+1,:], carr, rsmall, xg, yg, nr)

        k=k+1

        if k%1==0:
            # Clear the plot and redraw
            plt.clf()

            plotfield2(xg, yg, streamfield, streamname, xstore[k,:], ystore[k,:])

            plt.draw()
            plt.pause(0.001)
    plt.show()