from mesh_utils import mesh_output_writer
from mesh_utils import node_coord_displ
from mesh_utils import nodal_forces

#h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen = mesh_inputs()

h = 1.
l = 10.
m = 10
n = 100
option = 1.
E = 200000.
nu = 0.3
t = 1.
load_opt = 0
Mag = 0
R = 0
RadD = 0
V0 = 100.
filen = 'bob'

coord, displ, elnode, node, el = node_coord_displ(h, l, m, n, load_opt, R)

#forces, forcesb, DispMat, LoadMat = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

forces, forcesb = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

mesh_output_writer(filen, node, coord, el, option, elnode, E, nu, t)







#from numpy import array


#TEMPLATE = \
#"""------------------------------------------
#This log file contains:
#x1 \t x2 \t x3 \t f(x)
#------------------------------------------
#%s
#------------------------------------------"""

#LINE_ENTRY = "%.3f \t %.3f \t %.3f \t %.3f\n"

#def write_file(data):
    # initialize the string that contains the values from
    # the data array
#    body = ""
    # loop over each row vector in the data array (matrix)
#    for vector in data:
        # use the LINE_ENTRY string template to format the
        # values of the row vector into a string and join
        # this to body
#        body += LINE_ENTRY % tuple(vector)
        # remove the last newline charater
#        body = body.rstrip("\n")
        # populate the log file TEMPLATE with body
#        log_info = TEMPLATE % body
        # display the log file info to the screen
#        print log_info
        # open the data.txt file (write mode [w]) to write the
        # log_info
#        f_h = open("data.txt", "w")
        # write the log_info into the data.txt file
#        f_h.write(log_info)
        # close the data.txt file
#        f_h.close()
        
#def main():
    # create a data array for testing
#    data = array([[0.001, 0.011, 2.001, 10.123],
#                  [0.002, 0.021, 3.001, 21.103],
#     [0.003, 0.031, 4.001, 32.123],
#     [0.004, 0.041, 5.001, 43.123],
#     [0.005, 0.051, 6.001, 54.120],
#     [0.006, 0.061, 7.001, 65.123],
#     [0.007, 0.071, 8.001, 76.123],
#     [0.008, 0.081, 9.001, 87.023]], dtype=float)
     # call the function to write the data array to a file
#     write_file(data)
#     if __name__ == "__main__":
#         main()