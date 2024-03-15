"""Configuration file"""


############################
###        Domain        ###
############################

# Minimum latitude of domain
LAT_MIN = 61  # °N

# Latitude of sponge layer
LAT_SPONGE = 65  # °N


############################
###        Planet        ###
############################

# Planet radius
R = 66854e3  # m

# Gravity
G = 24.8  # m/s2

# Period
T = 9*3600 + 55*60 + 27  # s

# Rossby number
RO = 0.23

# Burger number 
BU = 10

############################
###  Initial conditions  ###
############################

# Steepness parameter
B = 1.5

# Vortex radius
R_M = 1000e3  # m

############################
###     Miscellaneous    ###
############################

# Courant number
CFL = 1.587 / 2

# Output timestep (in seconds or as a time string)
OUTPUT_TIMESTEP = '1 day'  # Note: it will be rounded to a multiple of dt

# Enable/disable logs
LOG = True

# Path to ffmpeg executable (None or empty string to keep matplotlib default value)
FFMPEG_PATH = R'C:\ffmpeg\bin\ffmpeg.exe'

# Output folder
OUTPUT_FOLDER = R'.\output'