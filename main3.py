# inheriting 2 classes with both containing a common attribute name may result in problems with MRO (multilpe resolution order) as ambiguity with which attribute youre reffering to exits 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import RegularPolygon
from matplotlib.patches import Circle
from math import pi
from math import sqrt 
from math import log
from math import isclose 
import random

fig, ax = plt.subplots()

users = []
basestations = []
monitored_basestations = []
monitored_users = []

HEXAGON_SIDES = 6
NETWORK_SIZE_X = 6000                               # Network length x-direction (m)
NETWORK_SIZE_Y = 6000                               # Network length y-direction (m)
CELL_RADIUS = 1000/sqrt(3)                          # bs-to-bs distance (Km)
NUMBER_OF_USERS = 5000                              # Number of Users to simulate
BS_TX_POWER = 40                                    # BS transmit power in Watts
BS_HEIGHT = 30                                      # BS Height in m
USER_HEIGHT = 1.6                                   # User Height in m
FREQUENCY = 2 * (10**9)                             # Operating frequency in GHz
C = 3 * (10**8)                                     # Speed of light
BANDWIDTH = 10 * (10**6)                            # Operating bandwidth in Hz
BOLTZMAN = 1.3806503 * (10**-23)                    # Boltzman's constant
TEMP = 297                                          # Operating temperature in Kelvin
THERMAL_NOISE = BOLTZMAN * TEMP * BANDWIDTH         # N = kTB 
NOISE_FIGURE = 8                                    # Noise Figure (NF) in dB 
NOISE_FACTOR = 10 ** (NOISE_FIGURE/10)              # Noise Factor (Fn) in decimal 
FRF = 3                                             # Frequency reuse Factor for FFR scheme
SCHEDULING_ROUNDS = 10                              # Number of scheduling rounds to simulate
INTERIOR_RADIUS_FACTOR = 0.6                        # Size of interior radius relative to cell radius (0-1 range) for FFR scheme


# dictionary with key-value pairs corresponding to basestation id: basestation of that bs
bsID_bs_dic = {}

class Coordinates:
    def __init__(self, x, y):
        self.x_coord = x
        self.y_coord = y

class User(Coordinates):
    def __init__(self, x, y):
        super().__init__(x, y)
        self.user_id = 0                        # default id 
        self.coords = (self.x_coord, self.y_coord)
        self.distance_to_all_bs = []            # this is a list of tuples, each containg the id of bs and distance to it 
        self.bs = 0                             # the bs class assigned to the user
        self.bs_coords = (0,0)                  # coordinates of the bs the user is connected to 
        self.first_tier_cells = []
        self.snr = 0
        self.sinr = 0
        self.inter_power = 0
        self.sig_power = 0
        self.capacity = 0                       # Shanon's capacity 
        self.rbs = []
        self.bw = 0
        self.tti = 0
        self.sinr_total = 0
        self.snr_total = 0
        self.interior = False
        self.exterior = False
    
    def plot_user(self):
        ax.plot(self.x_coord, self.y_coord, 'bo', markersize=0.5)
        # ax.annotate(self.user_id, (self.x_coord, self.y_coord))
        
    def distance_to_bs(self):
        distance_xy = calc_direct_distance(self.coords, self.bs_coords)
        return distance_xy

    def pl(self, distance, frequency):
        # return 30*log(distance, 10) + 30*log(frequency, 10) + 30*log(4*pi/C, 10)
        """ 3GPP path loss model (assumes 2GHz) """
        return 125.5 + 36.2*log(distance/1000, 10)
    
    def calc_snr(self):
        tx_power_dBW = 10*log(BS_TX_POWER, 10)
        distance_xy = calc_direct_distance(self.coords, self.bs_coords)
        distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
        received_power = tx_power_dBW - self.pl(distance_xyz, FREQUENCY)
        noise = 10*log(THERMAL_NOISE, 10) + NOISE_FIGURE
        snr_dB = received_power - noise
        self.snr = snr_dB
        return snr_dB
    
    def calc_sinr(self):
        tx_power_dBW = 10*log(BS_TX_POWER, 10)
        distance_xy = calc_direct_distance(self.coords, self.bs_coords)
        distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
        signal_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY) 
        interference = 0
        for bs_id in self.first_tier_cells:
            bs_coords = bsID_bs_dic[bs_id].coords
            distance_xy = calc_direct_distance(self.coords, bs_coords)
            distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
            received_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY)
            received_power = 10**(received_power_dBW/10)
            interference += received_power
        noise_dB = 10*log(THERMAL_NOISE, 10) + NOISE_FIGURE
        noise = 10**(noise_dB/10)
        noise_inter = noise + interference
        sinr_dB = signal_power_dBW - 10*log(noise_inter, 10)
        self.sinr = sinr_dB
        return sinr_dB
    
    def calc_sinr_ffr(self):
        tx_power_dBW = 10*log(BS_TX_POWER, 10)
        distance_xy = calc_direct_distance(self.coords, self.bs_coords)
        distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
        signal_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY) 
        noise_dB = 10*log(THERMAL_NOISE, 10) + NOISE_FIGURE
        noise = 10**(noise_dB/10)
        interference = 0

        if self.interior:
            for bs_id in self.first_tier_cells:
                bs_coords = bsID_bs_dic[bs_id].coords
                distance_xy = calc_direct_distance(self.coords, bs_coords)
                distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
                received_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY)
                received_power = 10**(received_power_dBW/10)
                interference += received_power
            noise_inter = noise + interference
            sinr_dB = signal_power_dBW - 10*log(noise_inter, 10)
            self.sinr = sinr_dB
            return sinr_dB
        
        if self.exterior:
            sinr_dB = signal_power_dBW - noise_dB
            self.sinr = sinr_dB
            return sinr_dB

    def calc_sinr_coop(self):
        tx_power_dBW = 10*log(BS_TX_POWER, 10)
        distance_xy = calc_direct_distance(self.coords, self.bs_coords)
        distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
        signal_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY) 
        interference = 0
        for bs_id in self.bs.first_tier_interferers:
            bs_coords = bsID_bs_dic[bs_id].coords
            distance_xy = calc_direct_distance(self.coords, bs_coords)
            distance_xyz = sqrt((distance_xy ** 2) + ((BS_HEIGHT - USER_HEIGHT) ** 2))
            received_power_dBW = tx_power_dBW - self.pl(distance_xyz, FREQUENCY)
            received_power = 10**(received_power_dBW/10)
            interference += received_power
        noise_dB = 10*log(THERMAL_NOISE, 10) + NOISE_FIGURE
        noise = 10**(noise_dB/10)
        noise_inter = noise + interference
        sinr_dB = signal_power_dBW - 10*log(noise_inter, 10)
        self.sinr = sinr_dB
        return sinr_dB
   
class Basestation(Coordinates):
    def __init__(self, x, y, radius):
        super().__init__(x, y)
        self.radius = radius
        self.basestation_id = 0             #default id
        self.coords = (self.x_coord, self.y_coord)
        self.first_tier_cells = []
        self.monitored = False 
        self.hexagon = 0
        self.users = []
        self.rbs = []
        self.interior_rbs = []
        self.first_tier_interferers = []
        self.scheduling_queue = []
        self.users_waiting = []
        self.tti = 0
        self.int_scheduling_queue = []
        self.ext_scheduling_queue = []

    def create_station(self):
        hex = RegularPolygon((self.x_coord, self.y_coord), HEXAGON_SIDES, self.radius, orientation=1.570796, facecolor='bisque', edgecolor='k', linestyle='solid', linewidth=2)
        self.hexagon = hex #update self.hexagon attribute 
        ax.add_patch(hex)
        ax.annotate(self.basestation_id, (self.x_coord, self.y_coord), size=10)

    def update_cell_colour(self):
        self.hexagon.set_alpha(0.2)
        ax.add_patch(self.hexagon)
    
    def color_cell(self):
        self.hexagon.set_color('red')
        ax.add_patch(self.hexagon)
    
# generator function that iterates between values 0 and 1 
def alternate():
    while True:
        yield 0
        yield 1

def calc_cells_x():
    cell_radius = CELL_RADIUS
    two_cell_radius = 2*CELL_RADIUS
    gen = alternate()
    i = gen.__next__()
    x_dist = 0
    num_cells_x = 0
    while True:
        x_iterator = [two_cell_radius, cell_radius]
        if (x_dist + x_iterator[i]) <= NETWORK_SIZE_X:            
            x_dist += x_iterator[i]
            #if the final cell doesn't completely fit within network then don't add this cell 
            if i == 1 and (x_dist + cell_radius/2 > NETWORK_SIZE_X):
                return num_cells_x
            num_cells_x += 1
            i = gen.__next__()
        else:             
            # print('num_cells_x: ', num_cells_x)
            return num_cells_x

def calc_cells_y():
    cell_height = sqrt(3)*CELL_RADIUS
    y_dist = 0
    num_cells_y = 0
    while True:
        if (y_dist + cell_height) <= NETWORK_SIZE_Y:
            y_dist += cell_height
            num_cells_y += 1
        else:
            # print('y_dist: ', y_dist)
            # print('num_cells_y: ', num_cells_y)
            return num_cells_y

def plot_network():
    ax.set_title('Network Simulation', fontsize=16)
    ax.set_xlabel('X distance (m)', fontsize=14)
    ax.set_ylabel('Y distance (m)', fontsize=14)
    
    for i_y in range(calc_cells_y()):
        gen = alternate()
        i = gen.__next__()      #0
        radius = CELL_RADIUS
        x = radius 
        y_iterator = [radius*sqrt(3)/2 + i_y*sqrt(3)*radius, radius*sqrt(3)+ i_y*sqrt(3)*radius]
        for i_x in range(calc_cells_x()):
            y = y_iterator[i]
            # print('i_y:{}, i_x:{}, x:{}, y:{}'.format(i_y, i_x, x, y))
            bs = Basestation(x, y, radius)
            bs.basestation_id = len(basestations)
            basestations.append(bs)
            bsID_bs_dic[bs.basestation_id] = bs
            bs.create_station()
            x += radius*3/2
            i = gen.__next__()

def find_monitored_cells():
    radius = CELL_RADIUS
    for bs in basestations:
        first_tier_stations = []
        x,y = bs.coords
        north_west = (x - 3/2*radius, y + sqrt(3)/2*radius)
        south_west = (x - 3/2*radius, y - sqrt(3)/2*radius)
        north = (x, y + sqrt(3)*radius)
        south = (x, y - sqrt(3)*radius)
        north_east = (x + 3/2*radius, y + sqrt(3)/2*radius)
        south_east = (x + 3/2*radius, y - sqrt(3)/2*radius)

        for i in basestations:
            if coords_is_close(i.coords, north_west):
                first_tier_stations.append(i.basestation_id)
            elif coords_is_close(i.coords, south_west):
                first_tier_stations.append(i.basestation_id)
            elif coords_is_close(i.coords, north):
                first_tier_stations.append(i.basestation_id)
            elif coords_is_close(i.coords, south):
                first_tier_stations.append(i.basestation_id)
            elif coords_is_close(i.coords, north_east):
                first_tier_stations.append(i.basestation_id)
            elif coords_is_close(i.coords, south_east):
                first_tier_stations.append(i.basestation_id)

        if len(first_tier_stations) == 6:
            monitored_basestations.append(bs) 
            bs.first_tier_cells = first_tier_stations
            bs.monitored = True
        
        # Update cell colour for non-moniotred cells
        else:
            bs.first_tier_cells = first_tier_stations
            bs.update_cell_colour()

def plot_users(num_users):
    current_num_users = 0
    while current_num_users < num_users:
        x = random.randint(0, NETWORK_SIZE_X)
        y = random.randint(0, NETWORK_SIZE_Y + int(sqrt(3)*CELL_RADIUS/2))
        point = (x,y)
        for bs in basestations:
            """ If the random point generated lies within an exisiting basestation, then a user is assigned to this point """
            if bs.hexagon.contains_point(ax.transData.transform(point)):
                user = User(x, y)
                user.user_id = len(users)
                user.bs = bs
                user.bs_coords = (bs.x_coord, bs.y_coord)
                user.first_tier_cells = bs.first_tier_cells
                if bs.monitored:
                    monitored_users.append(user)
                    # user.calc_snr()
                user.plot_user()
                users.append(user)
                bs.users.append(user)
                current_num_users += 1
                break

# Calculate Euclidean distance between 2 points
def calc_direct_distance(point1: tuple, point2: tuple):
    x1, y1 = point1
    x2, y2 = point2
    delta_x = x1 - x2
    delta_y = y1 - y2
    direct_distance = sqrt((delta_x * delta_x) + (delta_y * delta_y))
    return direct_distance   

#Returns True if 2 coordinates are equal to 4dp
def coords_is_close(coord1: tuple, coord2: tuple):
    if isclose(coord1[0], coord2[0], rel_tol=0.001) and isclose(coord1[1], coord2[1], rel_tol=0.001):
        return True 
    else:
        return False

def user_hotspot(num_users, radius):
    x = random.randint(int(2*CELL_RADIUS), int(NETWORK_SIZE_X - 2*CELL_RADIUS))
    y = random.randint(int(2*CELL_RADIUS), int(NETWORK_SIZE_Y - 2*CELL_RADIUS))
    bs_num = random.randint(0, len(basestations)-1)
    mon_bs = []
    for bs in monitored_basestations:
        mon_bs.append(bs.basestation_id)
    bs_num = random.choice(mon_bs)
    bs = basestations[bs_num]
    x,y = bs.coords
    coord = (x,y)
    hotspot = Circle(coord, radius)
    hotspot.set_facecolor('Red')
    hotspot.set_alpha(0.2)
    ax.add_patch(hotspot)
    current_users = 0
    while current_users < num_users:
        x1 = random.randint(int(x-radius), int(x+radius))
        y1 = random.randint(int(y-radius), int(y+radius))
        point = (x1, y1)
        if hotspot.contains_point(ax.transData.transform(point)):
            user = User(x1, y1)
            user.user_id = len(users)
            user.bs = bs
            user.bs_coords = (bs.x_coord, bs.y_coord)
            user.first_tier_cells = bs.first_tier_cells
            if bs.monitored:
                monitored_users.append(user)
            user.plot_user()
            users.append(user)
            bs.users.append(user)
            current_users += 1

def rr_scheduling(num_of_ttis):
    rb_size = 200*(10**3)
    rbs = int(BANDWIDTH/rb_size)
    scheduling_round = 0

    for bs in basestations:
        bs.scheduling_queue = [user for user in bs.users]
        bs.tti = num_of_ttis 

    while scheduling_round < num_of_ttis:
        for bs in basestations:
            for rb in range(1, rbs+1):
                active_user = bs.scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size
                bs.scheduling_queue.insert(0, active_user)

        calc_metrics()
        for bs in basestations:
            for user in bs.users:
                user.rbs = []
                user.bw = 0
        scheduling_round += 1
    calc_averages()

def bw_split(num_of_ttis, sinr_threshold):
    interior_users = 0
    for bs in monitored_basestations:
        bs.tti = num_of_ttis
        for user in bs.users:
            if user.calc_sinr() >= sinr_threshold:
                user.interior = True
                bs.int_scheduling_queue.append(user)
                interior_users += 1
            else:
                user.exterior = True
                bs.ext_scheduling_queue.append(user)
    
    print('Prop of interior users:', interior_users/len(monitored_users))
    return interior_users/len(monitored_users)

def ffr_scheduling(num_of_ttis, int_radius_factor):
    interior_radius = int(int_radius_factor*CELL_RADIUS)
    FRF = 3
    rb_size = 200*(10**3)

    interior_rbs = int((int_radius_factor**2)*(BANDWIDTH/rb_size))
    exterior_rbs = int((1-(int_radius_factor**2))*(BANDWIDTH/rb_size)/FRF)

    for bs in monitored_basestations:
        bs.tti = num_of_ttis 
        organise_users(bs, interior_radius)

    scheduling_round = 0
    while scheduling_round < num_of_ttis:
        for bs in monitored_basestations:
            for rb in range(1, interior_rbs+1):
                active_user = bs.int_scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size
                bs.int_scheduling_queue.insert(0, active_user)

            for rb in range(1, exterior_rbs+1):
                active_user = bs.ext_scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size
                bs.ext_scheduling_queue.insert(0, active_user)
        
        calc_metrics_ffr()
        for bs in monitored_basestations:
            for user in bs.users:
                user.rbs = []
                user.bw = 0
        scheduling_round += 1
    calc_averages()

def affr_scheduling(num_of_ttis):
    rb_size = 200*(10**3)
    
    split = bw_split(num_of_ttis, 9.5)
    interior_rbs = int(split*(BANDWIDTH/rb_size))
    exterior_rbs = int((1-split)*(BANDWIDTH/rb_size)/FRF)
 
    scheduling_round = 0
    while scheduling_round < num_of_ttis:
        for bs in monitored_basestations:
            for rb in range(1, interior_rbs+1):
                active_user = bs.int_scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size
                bs.int_scheduling_queue.insert(0, active_user)

            for rb in range(1, exterior_rbs+1):
                active_user = bs.ext_scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size
                bs.ext_scheduling_queue.insert(0, active_user)
        
        calc_metrics_ffr()
        for bs in monitored_basestations:
            for user in bs.users:
                user.rbs = []
                user.bw = 0
        scheduling_round += 1
    calc_averages()

def organise_users(bs, interior_radius):
    for user in bs.users:
        if user.distance_to_bs() < interior_radius:
            user.interior = True
            bs.int_scheduling_queue.append(user)
        
        elif user.distance_to_bs() >= interior_radius:
            user.exterior = True
            bs.ext_scheduling_queue.append(user)

def approach_1_1(num_of_ttis):
    rb_size = 200*(10**3)
    rbs = int(BANDWIDTH/rb_size)
    scheduling_round = 0

    # initial scheduling queue given by order of users sorted 
    for bs in basestations:
        bs.scheduling_queue = [user for user in bs.users]
        bs.tti = num_of_ttis

    rbs = [rb for rb in range(1, rbs+1)]
    i=0
    while scheduling_round < num_of_ttis:
        for bs in basestations:
            # clear bs.rbs from previous round
            bs.rbs = []
            
        for bs in (basestations):
            # find all rbs that have been assigned to first-tier cells of bs
            assigned_rbs = []
            for ft_bs in bs.first_tier_cells:
                assigned_rbs += bsID_bs_dic[ft_bs].rbs
          
            # find all rbs that can be assigned to the current bs (I.e., rbs not assigned to first-tier cells of bs)
            available_rbs = []
            for rb in rbs:
                if rb not in assigned_rbs:
                    available_rbs.append(rb)
  
            # append all available rbs as a bs attribute up to the number of users in bs
            bs.rbs = available_rbs[:len(bs.users)]
            # if not len(bs.rbs) == 0:
            #     bs.color_cell()
            for rb in bs.rbs:
                active_user = bs.scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size                  
                bs.scheduling_queue.insert(0, active_user)
        i += 1
        calc_metrics_coop()
        scheduling_round += 1
        for bs in basestations:
            for user in bs.users:
                user.rbs = []
                user.bw = 0
    calc_averages()

def approach_1_2(num_of_ttis):
    rb_size = 200*(10**3)
    rbs = int(BANDWIDTH/rb_size)
    scheduling_round = 0

    # initial scheduling queue given by order of users sorted 
    for bs in basestations:
        bs.scheduling_queue = [user for user in bs.users]
        bs.tti = num_of_ttis

    rbs = [rb for rb in range(1, rbs+1)]
    i=0
    while scheduling_round < num_of_ttis:
        for bs in basestations:
            # clear bs.rbs from previous round
            bs.rbs = []

        for bs in (basestations[i:] + basestations[:i]):        
            # find all rbs that have been assigned to first-tier cells of bs
            assigned_rbs = []
            for ft_bs in bs.first_tier_cells:
                assigned_rbs += bsID_bs_dic[ft_bs].rbs
            # print(f'Assigned RBs = {len(assigned_rbs)}')
            # find all rbs that can be assigned to the current bs (I.e., rbs not assigned to first-tier cells of bs)
            available_rbs = []
            for rb in rbs:
                if rb not in assigned_rbs:
                    available_rbs.append(rb)

            # append all available rbs as a bs attribute up to the number of users in bs
            bs.rbs = available_rbs[:len(bs.users)]
  
            for rb in bs.rbs:
                active_user = bs.scheduling_queue.pop()
                active_user.rbs.append(rb)
                active_user.bw += rb_size                  
                bs.scheduling_queue.insert(0, active_user)
        i += 1
        calc_metrics_coop()
        scheduling_round += 1
        for bs in basestations:
            for user in bs.users:
                user.rbs = []
                user.bw = 0
    calc_averages()
     
def calc_sinr():
    for user in monitored_users:
        user.calc_sinr()

def calc_metrics():
    for user in monitored_users:
        snr_dB = user.calc_snr()
        user.snr_total += snr_dB
        sinr_dB = user.calc_sinr()
        user.sinr_total += sinr_dB
        sinr_watts = 10**(sinr_dB/10)
        user.capacity += user.bw*log(1+sinr_watts, 2) #Shanon's capacity 

def calc_metrics_ffr():
    for user in monitored_users:
        snr_dB = user.calc_snr()
        user.snr_total += snr_dB
        sinr_dB = user.calc_sinr_ffr()
        user.sinr_total += sinr_dB
        sinr_watts = 10**(sinr_dB/10)
        user.capacity += user.bw*log(1+sinr_watts, 2) #Shanon's capacity 

def calc_metrics_coop():
    for user in monitored_users:
        snr_dB = user.calc_snr()
        user.snr_total += snr_dB
        sinr_dB = user.calc_sinr_coop()  
        # sinr_dB = user.calc_sinr_ffr() 
        user.sinr_total += sinr_dB
        sinr_watts = 10**(sinr_dB/10)
        user.capacity += user.bw*log(1+sinr_watts, 2) #Shanon's capacity 

def calc_averages():
    for user in monitored_users:
        user.capacity = user.capacity/user.bs.tti
        user.sinr = user.sinr_total/user.bs.tti
        user.snr = user.snr_total/user.bs.tti

def cdf(data):
    n = len(data)
    x = np.sort(data)
    y = np.arange(1, n + 1) / n 
    return x, y

def plot_snr_cdf():
    snr_list = []
    fig, ax = plt.subplots()
    for user in monitored_users:
        snr_list.append(user.snr)
    x_data,y_data = cdf(snr_list)
    ax.plot(x_data, y_data*100)
    ax.set_title('CDF of User SNR', fontsize=16)
    ax.set_xlabel("SNR (dB)", fontsize=15)
    ax.set_ylabel("CDF Percentage", fontsize=15)
    ax.grid(b=True)
    plt.show(block=False)

def plot_sinr_cdf():
    sinr_list = []
    fig, ax = plt.subplots()
    for user in monitored_users:
        sinr_list.append(user.sinr)
    x_data,y_data = cdf(sinr_list)
    ax.plot(x_data, y_data*100)
    ax.set_title('CDF of User SINR', fontsize=16)
    ax.set_xlabel("SINR (dB)", fontsize=15)
    ax.set_ylabel("CDF Percentage", fontsize=15)
    ax.grid(b=True)
    plt.show(block=False)
    
def plot_capacity_cdf():
    capacity_list = []
    fig, ax = plt.subplots()
    for user in monitored_users:
        capacity_list.append(user.capacity)
    capacity_mbps = [x / 10**6 for x in capacity_list]
    x_data,y_data = cdf(capacity_mbps)
    ax.plot(x_data, y_data*100)
    ax.set_title('CDF of User Capacity', fontsize=16)
    ax.set_xlabel("User Capacity (Mbps)", fontsize=15)
    ax.set_ylabel("CDF Percentage", fontsize=15)
    ax.grid(b=True)
    plt.show(block=False)

def plot_curves():
    plot_snr_cdf()
    plot_sinr_cdf() 
    plot_capacity_cdf()

plot_network()
find_monitored_cells()
plot_users(NUMBER_OF_USERS)
rr_scheduling(SCHEDULING_ROUNDS)
# approach_1_1(SCHEDULING_ROUNDS)
# approach_1_2(SCHEDULING_ROUNDS)
# ffr_scheduling(SCHEDULING_ROUNDS, INTERIOR_RADIUS_FACTOR)
# affr_scheduling(SCHEDULING_ROUNDS)
plot_curves()

""" Simulation results and variables """
def simulate_results():
    snr_total = 0
    sinr_total = 0
    capacity_total = 0
    not_scheduled = 0
    users_in_outage = 0
    for user in monitored_users:
        snr_total += user.snr
        sinr_total += user.sinr
        capacity_total += user.capacity
        if user.capacity < 0.1*(10**6):
            users_in_outage += 1
        if user.capacity == 0:
            not_scheduled += 1
    
    scheduled_users = len(monitored_users)-not_scheduled
    
    # Average values to 3 decimal places 
    snr_avg = "{:.3f}".format(snr_total / len(monitored_users))
    sinr_avg = "{:.3f}".format(sinr_total / len(monitored_users))
    # sinr_avg = "{:.3f}".format(sinr_total / scheduled_users)
    capacity_avg = "{:.3f}".format(capacity_total / (len(monitored_users) * (10**6)))
    capacity_avg_sched = "{:.3f}".format(capacity_total / ((scheduled_users) * (10**6)))
    outage_percentage = "{:.3f}".format((users_in_outage / len(monitored_users)) * 100)
    not_scheduled_percentage = "{:.3f}".format((not_scheduled / len(monitored_users)) * 100)
    avg_users_per_cell = "{:.3f}".format(len(monitored_users)/len(monitored_basestations))

    print(f'''
            Simulation Results:
            
            Number of monitored basestations: {len(monitored_basestations)}
            Number of monitored user: {len(monitored_users)}
            Average number of users per cell: {avg_users_per_cell}
            Average SNR: {snr_avg} dB
            Average SINR: {sinr_avg} dB
            Average Capacity of scheduled users: {capacity_avg_sched} Mbps
            User Outage: {outage_percentage}%
            Percentage of users not scheduled: {not_scheduled_percentage}%
            Total Network Capacity: {capacity_total / (10**6)} Mbps
        ''')
simulate_results()

plt.axis('equal')
plt.show() 


