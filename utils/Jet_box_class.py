import astropy.units as u
from astropy.coordinates import Angle
import math
import numpy as np
import matplotlib.pyplot as plt

def make_cluster_box(cluster_jet):
    
    ## calculate angle from subjects
    heights = []
    angles = []
    for jet in cluster_jet.jets:
        heights.append(jet.height)
        angles.append(jet.angle)
    jet_box_angle = angles[heights.index(max(heights))]
    jet_box_angle = np.mean(angles)
    
    jet_box = Jet_box([cluster_jet.Bx, cluster_jet.By], cluster_jet.Max_Height, cluster_jet.Width, jet_box_angle)
    return jet_box

def make_subject_box(subject_jet):
    subject_box = Jet_box(subject_jet.solar_start, subject_jet.solar_H, subject_jet.solar_W, subject_jet.angle)
    return subject_box

def plot_all_subject_boxes_and_average_box(jet_cluster, cluster_box, aia_map):
    
    fig = plt.figure()
    ax = fig.add_subplot(projection=aia_map)
    image = aia_map.plot(axes=ax)

    for subject in jet_cluster.jets:
        subject_box = Jet_box(subject.solar_start, subject.solar_H, subject.solar_W, subject.angle)
        for line in subject_box.lines_to_plot():
            ax.plot(line[0]*u.arcsec.to(u.deg), line[1]*u.arcsec.to(u.deg),
            color='white', linewidth=0.7, linestyle='dotted',
            transform=ax.get_transform("world"))

    for line in cluster_box.lines_to_plot():
        ax.plot(line[0]*u.arcsec.to(u.deg), line[1]*u.arcsec.to(u.deg),
            color='white',
            transform=ax.get_transform("world"))

    plt.show()

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
    
    def corners(self):
        x1 = self.base[0] - 0.5*self.width*np.sin(self.angle.radian)
        x2 = self.base[0] + 0.5*self.width*np.sin(self.angle.radian)
        y1 = self.base[1] + 0.5*self.width*np.cos(self.angle.radian)
        y2 = self.base[1] - 0.5*self.width*np.cos(self.angle.radian)
        dx = self.height*np.cos(self.angle.radian)
        dy = self.height*np.sin(self.angle.radian)
        dxp = self.width*np.sin(self.angle.radian)
        dyp = self.width*np.cos(self.angle.radian)
        x3 = x1+dx
        y3 = y1+dy
        x4 = x2+dx
        y4 = y2+dy
        return [[x1,y1],[x2,y2],[x3,y3],[x4,y4]]
    
    def lines_to_plot(self):
        corners = self.corners()
        line1x = np.array([corners[0][0].value, corners[2][0].value])
        line1y = np.array([corners[0][1].value, corners[2][1].value])
        line2x = np.array([corners[1][0].value, corners[3][0].value])
        line2y = np.array([corners[1][1].value, corners[3][1].value])
        line3x = np.array([corners[0][0].value, corners[1][0].value])
        line3y = np.array([corners[0][1].value, corners[1][1].value])
        line4x = np.array([corners[2][0].value, corners[3][0].value])
        line4y = np.array([corners[2][1].value, corners[3][1].value])
        return [[line1x, line1y], [line2x, line2y], [line3x, line3y], [line4x, line4y]]


