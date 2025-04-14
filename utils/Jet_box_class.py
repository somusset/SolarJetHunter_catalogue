import astropy.units as u
from astropy.coordinates import Angle
import math
import numpy as np
import matplotlib.pyplot as plt

"""
A few functions to create objects of the class Jet_box, or to plot the boxes on a sunpy map.
"""

def make_cluster_box(cluster_jet):
    
    ## calculate angle from subjects - instead of taking the average angle (makes little sense);
    ##                                  we keep the angle of the biggest box in the series of boxes - this is when the jet is most developed
    heights = []
    angles = []
    for jet in cluster_jet.jets:
        heights.append(jet.height)
        angles.append(jet.angle)
    jet_box_angle = angles[np.argmax(heights)]
    #jet_box_angle = np.mean(angles)
    
    jet_box = Jet_box([cluster_jet.Bx, cluster_jet.By], cluster_jet.Max_Height, cluster_jet.Width, jet_box_angle)
    return jet_box

def make_subject_box(subject_jet):
    subject_box = Jet_box(subject_jet.solar_start, subject_jet.solar_H, subject_jet.solar_W, subject_jet.angle)
    return subject_box

def plot_all_subject_boxes_and_average_box(jet_cluster, cluster_box, axes):
    for subject in jet_cluster.jets:
        subject_box = Jet_box(subject.solar_start, subject.solar_H, subject.solar_W, subject.angle)
        for line in subject_box.lines_to_plot():
            axes.plot(line[0]*u.arcsec.to(u.deg), line[1]*u.arcsec.to(u.deg),
            color='white', linewidth=0.7, linestyle='dotted',
            transform=axes.get_transform("world"))

    for line in cluster_box.lines_to_plot():
        axes.plot(line[0]*u.arcsec.to(u.deg), line[1]*u.arcsec.to(u.deg),
            color='white',
            transform=axes.get_transform("world"))
    plt.show()

"""
Jet_box class
"""

class Jet_box:
    # can be used for any box defined with solar coordinates
    
    def __init__(self, base, height, width, angle):
        self.base = base*u.arcsec
        self.height = height*u.arcsec
        self.width = width*u.arcsec
        self.angle = Angle(angle-math.pi/2, u.radian) 
        
    def area(self):
        return self.height*self.width
        
    def center(self):
        center_coordinates = [0.,0.]
        center_coordinates[0] = self.base[0] + self.height/2.*np.cos(self.angle.radian) 
        center_coordinates[1] = self.base[1] + self.height/2.*np.sin(self.angle.radian) 
        return center_coordinates
    
    def rotation_mat(self):
        # rotation matrix for the angle of the box
        rm = np.asarray([[np.cos(self.angle), -np.sin(self.angle)],
                    [np.sin(self.angle),  np.cos(self.angle)]])
        return rm

    def rotation_around_base(self, point):
        # point is a tuple or list or array or two quantities
        # this function returns the coordinates of the point after rotation around the base point of the box.
        # the rotation matrix is provided by another function
        point_x = point[0].to(u.arcsec).value
        point_y = point[1].to(u.arcsec).value
        unit = self.base[0].unit
        rotation_mat = self.rotation_mat()
        rot_inter = np.matmul(rotation_mat, (np.asarray([point_x,point_y]) - np.asarray([self.base[0].value, self.base[1].value])))
        rot_point = [self.base[0].value + rot_inter[0], self.base[1].value - rot_inter[1]]
        return rot_point*unit

    def corners(self, no_unit=False):
        # corners before rotation
        c1 = [self.base[0], self.base[1]-0.5*self.width]
        c2 = [self.base[0]+self.height, self.base[1]-0.5*self.width]
        c3 = [self.base[0]+self.height, self.base[1]+0.5*self.width]
        c4 = [self.base[0], self.base[1]+0.5*self.width]
        # rotation 
        rot_c1 = self.rotation_around_base(c1)
        rot_c2 = self.rotation_around_base(c2)
        rot_c3 = self.rotation_around_base(c3)
        rot_c4 = self.rotation_around_base(c4)
        if no_unit:
            result = [(rot_c1[0].value, rot_c1[1].value), (rot_c2[0].value, rot_c2[1].value), 
                      (rot_c3[0].value, rot_c3[1].value), (rot_c4[0].value, rot_c4[1].value)]
        else:
            result = [rot_c1, rot_c2, rot_c3, rot_c4]
        return result
    
    def lines_to_plot(self):
        corners = self.corners()
        line1x = np.array([corners[0][0].value, corners[1][0].value])
        line1y = np.array([corners[0][1].value, corners[1][1].value])
        line2x = np.array([corners[1][0].value, corners[2][0].value])
        line2y = np.array([corners[1][1].value, corners[2][1].value])
        line3x = np.array([corners[2][0].value, corners[3][0].value])
        line3y = np.array([corners[2][1].value, corners[3][1].value])
        line4x = np.array([corners[3][0].value, corners[0][0].value])
        line4y = np.array([corners[3][1].value, corners[0][1].value])
        return [[line1x, line1y], [line2x, line2y], [line3x, line3y], [line4x, line4y]]
