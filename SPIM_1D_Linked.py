#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:15:41 2023

@author: joshwolpert

1D model of river profile with linked tributaries. Simply run the script, and
respond to the prompted questions.
"""

# Import modules
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

# 1) Display Instructions with a function

# X) Call the main function main(), and encapsulate it in its own function
# and then just call the main function (which in turn calls the other functions)

# X) You'll want to write the output to a text file (See Chapter 7))

# X) If you write to a text file, have each timestep comprise a single line,
# then use:
#text_file=open("read_it.txt","r")
#print(text_file.readline())
# to read the entire line (or timestep in this case)

##############################################################################
# Define the Stream object
class Stream(object):
    """A Bedrock Stream"""
    
    # Constructor method
    def __init__(self, name):
        # Initialize stream object attributes
        self.name = name
        
        g = 0
        while g == 0:
            try:
                self.m = float(input("Area exponent (m): "))
                g = 1
            except ValueError:
                print("\nArea exponent must be a number.")
        
        g = 0
        while g == 0:
            try:
                self.n = float(input("Slope exponent(n): "))
                g = 1
            except ValueError:
                print("\nSlope exponent must be a number.")
                
        g = 0
        while g == 0:
            try:
                self.c = float(input("Hack coefficient (c): "))
                g = 1
            except ValueError:
                print("\nHack coefficient must be a number.")
                
        g = 0
        while g == 0:
            try:
                self.h = float(input("Hack exponent (h): "))
                g = 1
            except ValueError:
                print("\nHack exponent must be a number.")
                
        g = 0
        while g == 0:
            try:
                self.U = float(input("Initial rock uplift rate (Uo): "))
                g = 1
            except ValueError:
                print("\nUo must be a number.")
        
        g = 0
        while g == 0:
            try:
                self.K = float(input("Initial rock erodibility (Ko): "))
                g = 1
            except ValueError:
                print("\nKo must be a number.")
                
        g = 0
        while g == 0:
            try:
                self.l = float(input("Stream length (m): "))
                g = 1
            except ValueError:
                print("\nStream length must be a number.")
                
        g = 0
        while g == 0:
            try:
                self.dx = float(input("Node spacing (m): "))
                g = 1
            except ValueError:
                print("\nNode spacing must be a number.")
        
        # Check that stream length and node spacing are valid and compatible
        while np.mod(self.l, self.dx) or self.l==0 or self.dx==0:
            print("\nError: Node spacing must be a factor of the stream's length and neither node spacing nor stream length can equal 0.")
            self.l = float(input("\nPlease re-enter the stream's length (m): "))
            self.dx = float(input("\nPlease re-enter the stream's node spacing: "))
            
        self.nn = int(self.l/self.dx)+1
        
        g = 0
        while g ==0:
            try:
                self.Ai = float(input("Critical drainage area (Ai): "))
                g = 1
            except ValueError:
                print("\nCritical drainage area must be a number.")
            
        g = 0
        while g == 0:
            try:
                self.Ao = float(input("Reference drainage area for Chi calculation (Ao): "))
                g = 1
            except ValueError:
                print("\nReference drainage area must be a number.")

        self.x_scale = np.linspace(0,float(self.l),self.nn)

        self.A_scale = np.zeros(self.nn)
        self.e_scale = np.zeros(self.nn)
        self.ksn_scale = np.zeros(self.nn)
        
        # Set drainage area relationships with Hack's Law (Hack, 1957)
        for i in np.arange(len(self.A_scale)):
            self.A_scale[i] = self.Ai + (self.c*np.flip(self.x_scale)[i]**self.h)
        
        # Calculate chi of nodes along profile (Royden and Perron, 2013)
        self.integrand = np.power(np.divide(self.Ao,self.A_scale),(self.m/self.n))
        self.chi_scale = integrate.cumtrapz(self.integrand,dx=self.dx)
        self.chi_scale = np.insert(self.chi_scale,0,0) # Append 0 to downstream boundary of chi list
        
        # Calculate elevations of steady state profiles (Royden and Perron, 2013; Mitchell and Yanites, 2019)
        self.z_scale = np.multiply(self.chi_scale,0)
        self.zb = 0 # Downstream boundary condition
        
        for i in np.arange(len(self.z_scale)):
            if i == 0:
                self.z_scale[i] = self.zb
            else:
                self.z_scale[i] = self.zb + ((self.U/(self.K*(self.Ao**self.m)))**(1/self.n))*(self.chi_scale[i])
                
        # Set K_scale, which will be used to determine K at each position along a stream
        self.K_scale = np.zeros(self.nn)+self.K
        
        # Initialize list with times and uplift rates
        self.U_list = []


        ######################################################################
        # Enable user to change a Stream object's uplift rate
        answer = str(input("Would you like to add a change in uplift rate? (Yes/No): ")).upper()
        while answer != "YES" and answer != "NO":
            print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
            answer = str(input("Would you like to add a change in uplift rate? (Yes/No): ")).upper()
        
        while answer == "YES":
            g = 0
            while g == 0:
                try:
                    time = int(input("After how many years will a change in uplift rate occur?: "))
                    global tmax
                    global dt
                    if np.mod(time, dt) or time > tmax:
                        print("A change in uplift rate cannot occur after the model run has ended and must occur at a time step.")
                    else:
                        g = 1
                except ValueError:
                    print("Please enter a time (integer) that is a multiple of the model's time step (dx) and less than the model's run time (tmax)")
                
            g = 0  
            while g == 0:
                try:
                    rate = float(input("After " + str(time) + " years, what will the uplift rate change to?: "))
                    g = 1
                except ValueError:
                    print("The new uplift rate must be numeric. Please try again.")

            # Deposit times and rates of changes in U into list
            self.U_list.append([[time],[rate]])
            
            answer = str(input("Would you like to add another change in uplift rate? (Yes/No)")).upper()
            while answer != "YES" and answer != "NO":
                print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
                answer = str(input("Would you like to add a change in uplift rate? (Yes/No): ")).upper()
            
            # Sort changes in uplift rate by time
            np.sort(self.U_list, axis = 0)
                
                
        ######################################################################       
        # Enable user to define a Stream object's stratigraphy
        answer = str(input("Would you like to add a contact to the stream? (Yes/No): ")).upper()
        while answer != "YES" and answer != "NO":
            print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
            answer = str(input("Would you like to add a contact to the stream? (Yes/No): ")).upper()
            
        self.c_list = [] # Initialize list for starting depth at outlet and underlying K
        self.vert_c_list = [] # Initialize list for vertical contacts
        self.func_c_list = [] # Initialize list for function-defined contacts
        
        while answer == "YES": # User chooses to create a contact
            
            g = 0
            while g == 0:
                contact_type = str(input("Would you like to add a vertical contact, dipping planar contact, or a contact defined with a function? (Vertical/Dipping Planar/Define with Function): ")).upper()
                if contact_type == "VERTICAL" or contact_type == "DIPPING PLANAR" or contact_type == "DEFINE WITH FUNCTION":
                    g = 1
                else:
                    print("\nPlease try again and respond with either Vertical, Dipping Planar, or Define with Function.")
            
            # User chooses a vertical contact
            if contact_type == "VERTICAL":
                g = 0
                while g == 0:
                    try:
                        distance = float(input("At what upstream distance will this vertical contact outcrop?: "))
                        g = 1
                        if distance > self.l:
                            print("\nPlease try again. The contact's outcrop position needs to be a distance upstream from the stream's outlet and cannot exceed the stream's length.")
                            g = 0
                    except ValueError:
                        print("\nPlease try again. The contact's outcrop position needs to be a distance upstream from the stream's outlet and cannot exceed the stream's length.")

                g = 0
                while g == 0:
                    try:
                        contact_K = float(input("What is the erosional efficiency (K) upstream of the contact?: "))
                        g = 1
                    except ValueError:
                        print("\nK needs to be a real number value. Please try again.")
                
                # Append outcrop distances and upstream K of contacts to vert_c_list
                self.vert_c_list.append([[distance],[contact_K]])
                
                # Allow the user to add another contact
                answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
                while answer != "YES" and answer != "NO":
                    print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
                    answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
                
                
            if contact_type == "DIPPING PLANAR":
                g = 0
                while g == 0:
                    try:
                        elevation = float(input("At what depth/elevation will the contact intersect the outlet's along-stream position?: "))
                        g = 1
                    except ValueError:
                        print("\nPlease try again. The contact's outcrop depth/elevation must be a real number value.")
                        
                g = 0
                while g == 0:
                    try:
                        dip = float(input("What is the contact's dip angle in degrees? (Negative dips are towards the basin's divide, and Positive dips are towards the stream's outlet): "))
                        g = 1
                    except ValueError:
                        print("\nPlease try again. The contact's dip must be a real number value.")
                        
                g = 0
                while g == 0:
                    try:
                        contact_K = float(input("What is the erosional efficiency (K) underlying the contact?: "))
                        g = 1
                    except ValueError:
                        print("\nK needs to be a real number value. Please try again.")
                        
                # Append contact starting elevations at outlet, dip angles, and K underlying contacts to c_list
                self.c_list.append([[elevation],[dip],[contact_K]])
                
                # Allow the user to add another contact
                answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
                while answer != "YES" and answer != "NO":
                    print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
                    answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
                    
            
            if contact_type == "DEFINE WITH FUNCTION":
                g = 0
                while g == 0:
                    try:
                        strat_func = str(input("Enter a function defining the contact's dip, with f(0) equal to the contact's elevation/depth at the stream's outlet: "))
                        g = 1
                    except ValueError:
                        print("\nPlease try again. Function must include a variable (E.g., 2x**2 + 10")
                        
                g = 0
                while g == 0:
                    try:
                        contact_K = float(input("What is the erosional efficiency (K) underlying the contact?: "))
                        g = 1
                    except ValueError:
                        print("\nK needs to be a real number value. Please try again.")
                
                # Append function defining contact, and K underlying contacts to func_c_list
                self.func_c_list.append([[strat_func],[contact_K]])
                        
                # Allow the user to add another contact
                answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
                while answer != "YES" and answer != "NO":
                    print("\nPlease respond with either 'Yes' or 'No' and press enter to continue.")
                    answer = str(input("Would you like to add another contact to the stream? (Yes/No): ")).upper()
        
   
    # Plot stream
    def plot_profile(self):
        if ii == 0:
            plt.plot(self.x_scale, self.z_scale, 'b',label = 'Trunk')
            ax1.grid(which='major', axis='both')
        else:
            plt.plot(self.x_scale, self.z_scale, 'c',label = 'Tributary')
        ax1.set(xlabel='Distance from Outlet (m)', ylabel='Elevation (m)',
                title=str(t)+' years after trunk loses 25% drainage area')
        ax1.legend(framealpha = 1)
        
    def plot_erosion(self):
        if ii == 0:
            plt.plot(self.x_scale[1:], (self.e_scale*(1e6/dt))[1:], 'b',label = 'Trunk')
            ax1.grid(which='major', axis='both')
        else:
            plt.plot(self.x_scale[1:], (self.e_scale*(1e6/dt))[1:], 'c',label = 'Tributary')
        ax1.set(xlabel='Distance from Outlet (m)', ylabel='Erosion Rate (m)',
                title=str(t)+' years after trunk loses 25% drainage area')
        ax1.legend(framealpha = 1)
        
    # Plot contact
    def plot_contact(self,ii):
        # Check if trunk stream
        if ii == 0:
            for i in range(len(self.dipping_c_z)):
                cut_x = []
                cut_z = []
                for j in range(len(self.x_scale)):
                    if self.z_scale[j] > self.dipping_c_z[i][j]:
                        cut_x.append(self.x_scale[j])
                        cut_z.append(self.dipping_c_z[i][j])
                plt.plot(cut_x, cut_z,'r')
                ax1.grid()
            
            for j in range(len(self.vert_c_list)):
                ### Need to make sure plotting can happen.... location of vert contact must be on x_scale
                res = [idx for idx, val in enumerate(self.x_scale) if val == self.vert_c_list[j][0][0]]
                plt.plot([self.vert_c_list[j][0][0],self.vert_c_list[j][0][0]],[0, self.z_scale[res]],'r')
                ax1.grid()
        
        # Otherwise, it's a tributary. Don't plot beneath trunk ###### Working on this#######
        else:
            for i in range(len(self.dipping_c_z)):
                cut_x = []
                cut_z = []
                for j in range(len(self.x_scale)):
                    self.d = np.add(j,self.outlet_location/streams[0].dx)
                    if self.z_scale[j] > self.dipping_c_z[i][j] and self.dipping_c_z[i][j] > streams[0].z_scale[int(self.d)]:
                        cut_x.append(self.x_scale[j])
                        cut_z.append(self.dipping_c_z[i][j])
                plt.plot(cut_x, cut_z,'r')
                ax1.grid()
            
            for j in range(len(self.vert_c_list)):
                res = [idx for idx, val in enumerate(self.x_scale) if val == self.vert_c_list[j][0][0]]
                plt.plot([self.vert_c_list[j][0][0],self.vert_c_list[j][0][0]],[0, self.z_scale[res]],'r')
                ax1.grid()
    
    # Erode stream using explicit finite difference scheme
    def erode_profile(self):
        for i in range(len(self.z_scale)):
            if i != 0:
                self.e_scale[i] = self.K_scale[i]*(self.A_scale[i]**self.m)*(((self.z_scale[i]-self.z_scale[i-1])/self.dx)**self.n)*dt
        self.z_scale = self.z_scale - self.e_scale

                
    # Uplift a trunk stream
    def uplift_profile(self):
        self.z_scale[1:-1] = self.z_scale[1:-1]+(self.U*dt)
        
    
    # Check Courant condition for stability. Reduce dt by half if a stream becomes unstable
    def check_stability(self):
        a = np.multiply(np.multiply(-self.K_scale,np.power(self.A_scale,self.m)),np.power(np.insert(self.z_scale[1:len(self.z_scale)]-self.z_scale[0:len(self.z_scale)-1],0,0), self.n))
        global dt
        if np.max(np.divide(abs(a)*dt,self.dx)) > 0.9:
            dt /= 2
                    
                    
    # Check if U should change and make the change
    def check_U(self):
        if len(self.U_list) > 0:
            global t
            if self.U_list[0][0][0] == t:
                self.U = self.U_list[0][1][0]
                del self.U_list[0]
                np.sort(self.U_list, axis = 0)
            
    
    # Establish initial conditions for dipping contacts
    def contact_initial_conditions(self):
        # Populate lists for dipping contact elevations and K values
        self.dipping_c_z = []
        self.dipping_c_K = []
                
        # Iterate through planar contacts
        for i in range(len(self.c_list)):
            angle = np.tan(self.c_list[i][1][0]*np.pi/180)
            self.dipping_c_z.append(self.c_list[i][0][0]+self.z_scale[0]+((self.x_scale-self.x_scale[0])*angle))
            self.dipping_c_K.append([self.c_list[i][2]])

        # Iterate through function-defined contacts
        ####### I probably need to add a try/except statement to handle the eval function #############
        for i in range(len(self.func_c_list)):
            self.dipping_c_z.append(self.z_scale[0]+eval(self.func_c_list[i][0][0].replace('x','self.x_scale-self.x_scale[0]')))
            self.dipping_c_K.append([self.func_c_list[i][1]])
            
    # Uplift contacts
    def uplift_contacts(self):
        for i in range(len(self.dipping_c_z)):
            self.dipping_c_z[i] = self.dipping_c_z[i]+(self.U*dt)
            
    # Determine K_scale
    def set_K_scale(self):
        if len(self.vert_c_list) > 0:
            # First, iterate through vertical contacts and adjust K_scale
            self.vert_c_list = np.sort(self.vert_c_list, axis = 0)
            for i in range(len(self.vert_c_list)):
                res = [idx for idx, val in enumerate(self.x_scale) if val > self.vert_c_list[i][0][0]] 
                self.K_scale[res] = self.vert_c_list[i][1][0]
        
        # Assess dipping contacts
        if len(self.dipping_c_z) > 0:
            z_diffs = []
            for i in range(len(self.dipping_c_z)):
                z_diffs.append(self.dipping_c_z[i] - self.z_scale) # Subtract profile elevations from each contact's elevations
            z_diffs = np.matrix(z_diffs) # convert to a numpy matrix
            z_diffs[z_diffs < 0] = np.nan # Make negative contact depths nan
            
            for i in range(len(self.K_scale)):
                try:
                    d = np.nanargmin(z_diffs[:,i],0)
                    d = int(d[0,0])
                    self.K_scale[i] = self.dipping_c_K[d][0][0]
                except:
                    pass
                
    def calculate_ksn(self):
        for i in range(len(self.z_scale)):
            if i != 0:
                self.ksn_scale[i] = ((self.z_scale[i]-self.z_scale[i-1])/self.dx)/((self.A_scale[i])**((-1*self.m)/self.n))
        
##############################################################################
# Set up stream network

# Assign time info. Ensure dt is compatible with tmax
print("\nLet's begin by establishing temporal parameters.")

f = 0
while f == 0:
    g = 0
    while g == 0:
        try:
            dt = int(input("What is the model's time step duration (dt)?: "))
            g = 1
        except ValueError:
            print("\ndt needs to be an integer. Please try again.")
 
    g = 0
    while g == 0:
        try:
            tmax = int(input("What is the model's total run time (years)?: "))
            g = 1
        except ValueError:
            print("\ntmax needs to be an integer. Please try again.")
            
    if np.mod(tmax,dt):
        print("\ndt must be a factor of tmax. Please re-enter values.")
        f = 0
    else:
        f = 1 

# Store original dt for use during stability check
original_dt = dt

# Instantiate stream object to make the trunk stream
print("\nGreat! Now we'll construct a stream network! We'll start with the trunk stream.")
trunk = Stream("trunk")

g = 0
while g == 0:
    try:
        n_tribs = int(input("How many tributaries would you like to add?: "))
        g = 1
    except ValueError:
        print("Please provide an integer value.")

# Store stream objects in dictionary
streams = {}
streams[0] = trunk # 0th position is trunk stream

# Initialize list to store tributary confluence indices
outlet_idx = np.zeros(n_tribs+1)

tt = 1

# Create tributary network
while tt <= n_tribs:
    
    print("\nGreat! Now let's make tributary " + str(tt) + ".")
    
    g = 0
    while g == 0:
        try:
            outlet_location = int(input("Where is the confluence of this tributary and the trunk stream? Enter a distance upstream from the trunk's outlet: "))
            g = 1
            if outlet_location > np.max(streams[0].l):
                print("\nThe tributary's confluence location with the trunk stream cannot exceed the trunk stream's total length. Please try again.")
                g = 0
        except ValueError:
            print("The upstream confluence distance must be an integer. Please try another value.")
            
    
    # Instantiate tributary stream objects
    streams[tt] = Stream(tt)
    
    # Store outlet location
    streams[tt].outlet_location = outlet_location
    
    # Translate x_scale to account for upstream distance of outlet
    streams[tt].x_scale = streams[tt].x_scale + outlet_location
    
    # Find index of tributary confluence on trunk
    outlet_idx[tt] = streams[0].x_scale.tolist().index(outlet_location)
    
    # Translate z_scale to 'link' tributary with trunk stream
    streams[tt].z_scale = streams[tt].z_scale + streams[0].z_scale[int(outlet_idx[tt])]
    
    # Add upstream drainage area from tributary to trunk stream
    streams[0].A_scale[0:int(outlet_idx[tt])+1] = streams[0].A_scale[0:int(outlet_idx[tt])+1] + max(streams[tt].A_scale)
    
    tt += 1

# Construct initial contact elevations and set K_scale
for i in range(len(streams)):
    streams[i].contact_initial_conditions()
    streams[i].set_K_scale()


##############################################################################
# Main loop

# Initialize variables used in main loop
t = 0
r = 0
q = 1
A_change = 0

while t < tmax:
    
    # Return to original dt
    dt = original_dt
    
    # Check all streams for stability. If one stream is unstable, dt is reduced for all.
    for k in range(len(streams)):
        streams[k].check_stability()
    
    # Uplift and erode trunk stream
    streams[0].check_U()
    streams[0].uplift_profile()
    streams[0].erode_profile()
    
    ####### Edits to Simulate Change in Drainage Area Along Trunk Stream #######
    if t > 5000000:
        # Reset drainage area of trunk stream to original value
        streams[0].A_scale -= A_change
        # Calculate ksn once drainage area has been reset
        for i in range(len(streams)):
            streams[i].calculate_ksn()
        # Milankovitch timescale periods (~40 kyr). The sine of change in drainage area matters.
        # 40 kyr - lambda = 1.5708e-4,
        # 100 kyr - lambda = 6.2832e-5,
        A_change = -500000*np.sin(1.5708e-4*(t-5000000-10000))-500000 # -10000, -25000
        streams[0].A_scale += A_change
    ############################################################
    
    # n=1 - K = 1e-6
    # n=3 - K = 2.5e-11
    # n=0.8 - K = 2.89e-6
    
    # Link tributaries to trunk by making the 0th position of tributary z_scales corresponding trunk elevations
    for k in range(len(streams)):
        if k != 0:
            
            streams[k].z_scale[0] = streams[0].z_scale[int(outlet_idx[k])]
            
            # Uplift and erode tributaries.
            streams[k].check_U()
            streams[k].uplift_profile()
            streams[k].erode_profile()
            
            # Add the outlet node erosion rate to each tributary
            streams[k].e_scale[0] = streams[0].e_scale[int(outlet_idx[k])]
                    
    # Uplift contacts and set K_scale
    for i in range(len(streams)):
        streams[i].uplift_contacts()
        streams[i].set_K_scale()
        
    t += dt
    
    # Specify when and what to save
    if t>4999999:
        if t % 5000 == 0:
            for ii in range(len(streams)):
                np.savetxt("/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/Run10/st" + str(ii) + "_e_" + str(t) + ".csv",streams[ii].e_scale,delimiter=",")
                np.savetxt("/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/Run10/st" + str(ii) + "_z_" + str(t) + ".csv",streams[ii].z_scale,delimiter=",")
                np.savetxt("/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/Run10/st" + str(ii) + "_ksn_" + str(t) + ".csv",streams[ii].ksn_scale,delimiter=",")
                
# Plot model data
fig_1, ax1 = plt.subplots(1,1)
for ii in range(len(streams)):
     streams[ii].plot_profile()
     streams[ii].plot_contact(ii)
     streams[ii].plot_erosion
    
print(t)
