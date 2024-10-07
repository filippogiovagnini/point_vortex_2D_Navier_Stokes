import matplotlib.pyplot as plt
import numpy as np

# Compute the velocity at a point (x, y) due to a set of point vortices
def Velocity_at_point(xarr, yarr, carr, rsmall, x, y):
    dy = 0.05; dx = 0.05
    u_vel = -(streamfunction(xarr, yarr, carr, rsmall, x, y+dy) - streamfunction(xarr, yarr, carr, rsmall, x, y-dy))/2*dy
    v_vel =  (streamfunction(xarr, yarr, carr, rsmall, x+dx, y) - streamfunction(xarr, yarr, carr, rsmall, x-dx, y))/2*dx
    return u_vel, v_vel

# Update the position of the point vortices using the Forward Euler method
def FE(xarr, yarr, carr, rsmall):
    dt= 1
    for i in range(0,len(xarr)):
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        xarr[i] = xarr[i] + dt*u_vel
        yarr[i] = yarr[i] + dt*v_vel
    return xarr, yarr


def Leapfrog(xarr, yarr,xarr_old,yarr_old, carr, rsmall):
    dt = 1
    for i in range(0,len(xarr)):
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        xarr[i] = xarr_old[i] + 2*dt*u_vel
        yarr[i] = yarr_old[i] + 2*dt*v_vel
    return xarr, yarr
    
def Wiki(xarr, yarr, carr, rsmall, x, y):
    print("ncf, and probably wrong.")
    dt = 1
    w0 = -2**1/3/(2 - 2**1/3)
    w1 = 1/(2 - 2**1/3)
    c_coef = [w1/2,(w0+w1)/2,(w0+w1)/2,w1/2]
    d_coef = [w1,w0,w1,0]
    for i in range(0,len(xarr)):
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        xarr[i] = xarr[i] + c_coef[0]*dt*u_vel
        u_vel, v_vel = Velocity_at_point(xarr, yarr, carr, rsmall, xarr[i], yarr[i])
        yarr[i] = yarr[i] + d_coef[0]*dt*v_vel
            
    return xarr, yarr



def RK2(xarr, yarr, carr, rsmall, x, y):
    dt= 1
    xarr0, yarr0 = FE(xarr, yarr, carr, rsmall, x, y)
    xarr1, yarr1 = FE(xarr0, yarr0, carr, rsmall, x, y)
    xarr = 0.5*(xarr1 + xarr0)
    yarr = 0.5*(yarr1 + yarr0)
    return xarr, yarr

    
def Update_stream_field(xarr, yarr, carr, rsmall, xg, yg, nr):

    for j in range(ng):
        for i in range(ng):
            streamx = streamfunction(xarr, yarr, carr, rsmall, xg[i], yg[j])
            streamfield[i, j] = streamx
    return streamfield
    

# Compute the streamfunction at a point (x, y) due to a set of point vortices
def streamfunction(xarr, yarr, carr, rsmall, x, y):
    nv = len(carr)
    streamx = 0
    for i in range(nv):
        streamx = streamx + carr[i]*greenfunc(x, y, xarr[i], yarr[i], rsmall)
    return streamx


# def greenfunc(x, y, xi, yi, rsmall):
#     '''
#     Calculate the value of the Green function at position (x, y) on plane
#     given the centre of a point vortex is at (xi, yi).
#     Note that in special situation that (x, y)=(xi, yi) the Green function
#     cannot be evaluated this way. So for simplicity the distance separating
#     the points is set to a small positive value, rsmall.
    
#     Input: x, y, xi, yi, rsmall
#     Output: greenx - value of the Green function
#     '''
#     #
#     # Calculate the distance-squared separating (x, y) and (xi, yi)
#     #
#     rsq = (x-xi)**2 + (y-yi)**2
#     rsqpos = max([rsq, rsmall])
#     greenx = (1.0/4*np.pi)*np.log(rsqpos)
#     #greenx = (1.0/4*np.pi)*1/rsqpos**0.5
#     return greenx

def greenfunc(x, y, xi, yi, rsmall):
    rsq = (x-xi)**2 + (y-yi)**2
    rsqpos = max([rsq, rsmall])
    a = 0.0005
    p = 4
    greenx = (1.0/4*np.pi)*1/p*np.log(a*rsqpos**p + 1) + (1.0/4*np.pi)*np.log(rsqpos)
    return greenx

    
def plotfield(xg, yg, field, fieldname):
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
    nlevs = 120  # number of contour levels to use
    #
    # Note that plt.contourf plots an array as transposed due to convention on 
    # first index as columns and second as rows. So to get a
    # plot of field in expected orientation need to transpose the array (using .T).
    plt.contourf(xg, yg, field.T, nlevs)
    plt.show()

    return


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
    nlevs = 20  # number of contour levels to use
    #
    # Note that plt.contourf plots an array as transposed due to convention on 
    # first index as columns and second as rows. So to get a
    # plot of field in expected orientation need to transpose the array (using .T).
    plt.contourf(xg, yg, field.T, nlevs)
    plt.scatter(xarr, yarr, color='black')
    return
    


if __name__ == '__main__':

    #print('value of velocity at x, y = ',velpt)

    # Loop over a grid of points on the (x, y) plane to find the
    # streamfunction field and then plot it.
    # nr*rscal is the half width of the square grid domain chosen.
    #

    rsmall = 1.e-10
    nr = 10
    drscal = 0.2
    ng = 2*nr+1
    xg = drscal*(np.arange(ng)-nr)
    yg = np.copy(xg)
    streamfield = np.zeros([ng, ng])
    streamname = 'streamfunction'

    nv=11
    theta = np.linspace(0,2*np.pi,nv)
    xarr = np.cos(theta)
    yarr = np.sin(theta)
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

    while (k<nt-1):
        
        #xstore[k+1,:], ystore[k+1,:]= RK2(xstore[k,:], ystore[k,:], carr, rsmall, x, y)

        xstore[k+1,:], ystore[k+1,:]= Leapfrog(xstore[k,:], ystore[k,:],xstore[k-1,:], ystore[k-1,:], carr, rsmall)
        streamfield = Update_stream_field(xstore[k+1,:],ystore[k+1,:], carr, rsmall, xg, yg, nr)

        k=k+1

        if k%1==0:
            plt.clf()

            plotfield2(xg, yg, streamfield, streamname, xstore[k,:], ystore[k,:])
            plt.draw()
            plt.pause(0.001)
    plt.show()
    
    # while (k<nt-1):
#         k = k+1
#         xarr, yarr= FE(xarr, yarr, carr, rsmall, x, y)
#
#
#         streamfield = Update_stream_field(xarr, yarr, carr, rsmall, x, y,xg,yg)
#
#         xstore[k,:] = xarr
#         ystore[k,:] = yarr
#         plt.clf()
#         if (nt%100==0):
#             plotfield2(xg, yg, streamfield, streamname,xarr,yarr)
#             plt.draw()
#             plt.pause(0.001)
#     plt.show()
#
#     plt.clf()
#     plt.scatter(xstore,ystore)
#     plt.show()
    
    #u_vel, v_vel = Compute_Velocity(streamfield,drscal)


def VelocityUV(xarr,yarr,carr):
    nv = len(carr)
    u_vel = torch.zeros(nv)
    v_vel = torch.zeros(nv)
    numerator1 = torch.zeros(nv)
    numerator2 = torch.zeros(nv)
    denominator1 = torch.zeros(nv)
   
    alpha = 0.01
    delta  = 1e-3 #0.0001
   
    ## EULER
    for i in range(nv):
        denominator1 = (xarr-xarr[i])**2 + (yarr-yarr[i])**2
        u_vel += -carr[i]*1/2*np.pi*(yarr-yarr[i])*1 / denominator1
        v_vel +=  carr[i]*1/2*np.pi*(xarr-xarr[i])*1 / denominator1
 
   
 
    return u_vel, v_vel

def RK4(xarr, yarr, carr, dt):
    """deterministic rk4"""
 
    u_vel1, v_vel1 = VelocityUV(xarr,yarr,carr)
    xarr1 = xarr + dt*u_vel1/2
    yarr1 = yarr + dt*v_vel1/2
   
    u_vel2, v_vel2 = VelocityUV(xarr1,yarr1,carr)
    xarr2 = xarr + dt*u_vel2/2
    yarr2 = yarr + dt*v_vel2/2
    
    u_vel3, v_vel3 = VelocityUV(xarr2,yarr2,carr)
    xarr3 = xarr + dt*u_vel3
    yarr3 = yarr + dt*v_vel3
   
    u_vel4, v_vel4 = VelocityUV(xarr3,yarr3,carr)
    xarr4 = xarr + dt/6*(u_vel1+ 2*u_vel2  + 2*u_vel3 + u_vel4)
    yarr4 = yarr + dt/6*(v_vel1+ 2*v_vel2  + 2*v_vel3 + v_vel4)
   
    return xarr4, yarr4