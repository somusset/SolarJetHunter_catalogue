import json
from shapely.geometry import Polygon
#import astropy.units as u
#from astropy.coordinates import Angle
#import math
import numpy as np

def json_import_list(input_file):
    '''
        import a list of JetCluster objects from the input_file file.
        Inputs
            ------
            input_file : string
                path or filename to the json file with JetCluster objects
        Outputs
            ------
            clusters : list
                list of JetCluster objects
    '''
    with open(input_file, 'r') as file:
        lists = json.load(file)

    clusters = []

    for k in range(len(lists)):
        json_obj = lists[k]
        jets_subjson = json_obj['jets']

        jets_list = []

        for J in jets_subjson:
            subject = J['subject']
            best_start = np.array([J['start'][i] for i in ['x', 'y']])
            best_end = np.array([J['end'][i] for i in ['x', 'y']])
            jet_params = np.array([J['cluster_values'][i]
                                  for i in ['x', 'y', 'w', 'h', 'a']])
            jeti = Polygon(get_box_edges(*jet_params))
            jet_obj = Jet(subject, best_start, best_end, jeti, jet_params)
            jet_obj.time = np.datetime64(J['time'])
            jet_obj.sigma = J['sigma']
            
            if 'solar_cluster_values' in J:
                jet_obj.solar_cluster_values = np.array([J['solar_cluster_values'][i]
                                  for i in ['x', 'y', 'w', 'h', 'a']])
            jet_obj.solar_H = J['solar_H']
            jet_obj.solar_H_sig = np.array(
                [J['solar_H_sig'][i] for i in ['upper', 'lower']])
            jet_obj.solar_W = J['solar_W']
            jet_obj.solar_start = np.array(
                [J['solar_start'][i] for i in ['x', 'y']])
            jet_obj.solar_end = np.array(
                [J['solar_end'][i] for i in ['x', 'y']])
            jets_list.append(jet_obj)

        jets_list = np.asarray(jets_list)

        cluster_obj = JetCluster(jets_list)
        cluster_obj.ID = json_obj['id']
        cluster_obj.SOL = json_obj['SOL']
        cluster_obj.Duration = json_obj['duration']
        cluster_obj.obs_time = np.datetime64(json_obj['obs_time'])

        cluster_obj.Bx = json_obj['Bx']['mean']
        cluster_obj.std_Bx = json_obj['Bx']['std']

        cluster_obj.By = json_obj['By']['mean']
        cluster_obj.std_By = json_obj['By']['std']

        cluster_obj.Lat = json_obj['lat']
        cluster_obj.Lon = json_obj['lon']

        cluster_obj.Max_Height = json_obj['max_height']['mean']
        try:
            cluster_obj.std_maxH = np.array(
                [json_obj['max_height'][i] for i in ['std_upper', 'std_lower']])
        except Exception as e:
            print(e)
            cluster_obj.std_maxH = np.array([np.nan, np.nan])

        cluster_obj.Width = json_obj['width']['mean']
        cluster_obj.std_W = json_obj['width']['std']
        cluster_obj.Height = json_obj['height']['mean']
        cluster_obj.std_H = json_obj['height']['std']

        cluster_obj.sigma = json_obj['sigma']

        if 'velocity' in json_obj:
            cluster_obj.Velocity = json_obj['velocity']
        else:
            cluster_obj.Velocity = np.nan

        if 'flag' in json_obj:
            cluster_obj.flag = json_obj['flag']

        clusters.append(cluster_obj)

    clusters = np.asarray(clusters)

    print(f'The {len(clusters)} JetCluster objects are imported from {input_file}.')

    return clusters

def get_box_edges(x, y, w, h, a):
    '''
        Return the corners of the box given one corner, width, height
        and angle

        Inputs
        ------
        x : float
            Box left bottom edge x-coordinate
        y : float
            Box left bottom edge y-coordinate
        w : float
            Box width
        h : float
            Box height
        a : flat
            Rotation angle

        Outputs
        --------
        corners : numpy.ndarray
            Length 4 array with coordinates of the box edges
    '''
    cx = (2*x+w)/2
    cy = (2*y+h)/2
    centre = np.array([cx, cy])
    original_points = np.array(
        [
            [cx - 0.5 * w, cy - 0.5 * h], # This would be the box if theta = 0
            [cx + 0.5 * w, cy - 0.5 * h],
            [cx + 0.5 * w, cy + 0.5 * h],
            [cx - 0.5 * w, cy + 0.5 * h],
            # repeat the first point to close the loop
            [cx - 0.5 * w, cy - 0.5 * h]
        ]
    )
    rotation = np.array([[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]])
    corners = np.matmul(original_points - centre, rotation) + centre
    return corners

class Jet:
    '''
        Oject to hold the data associated with a single jet.
        Contains the start/end positions and associated extracts,
        and the box (as a `shapely.Polygon` object) and corresponding
        extracts
    '''

    def __init__(self, subject, start, end, box, cluster_values):
        self.subject = subject
        self.start = start
        self.end = end
        self.box = box

        self.cluster_values = cluster_values

        self.box_extracts = {'x': [], 'y': [], 'w': [], 'h': [], 'a': []}
        self.start_extracts = {'x': [], 'y': []}
        self.end_extracts = {'x': [], 'y': []}

        self.autorotate()
        
    def adding_new_attr(self, name_attr,value_attr):
        '''
            Add an additional attribute of value value_attr and name name_attr to the jet object 
        '''
        setattr(self, name_attr, value_attr)

    def get_extract_starts(self):
        '''
            Get the extract coordinates associated with the
            starting base points

            Outputs
            -------
            coords : numpy.ndarray
                coordinates of the base points from the extracts
                for the start frame
        '''
        x_s = self.start_extracts['x']
        y_s = self.start_extracts['y']

        return np.transpose([x_s, y_s])

    def get_extract_ends(self):
        '''
            Get the extract coordinates associated with the
            final base points

            Outputs
            -------
            coords : numpy.ndarray
                coordinates of the base points from the extracts
                for the final frame
        '''
        x_e = self.end_extracts['x']
        y_e = self.end_extracts['y']

        return np.transpose([x_e, y_e])

    def get_extract_boxes(self):
        '''
            Get the extract shapely Polygons corresponding
            to the boxes

            Outputs
            -------
            boxes : list
                List of `shapely.Polygon` objects corresponding
                to individual boxes in the extracts
        '''
        boxes = []
        for i in range(len(self.box_extracts['x'])):
            x = self.box_extracts['x'][i]
            y = self.box_extracts['y'][i]
            w = self.box_extracts['w'][i]
            h = self.box_extracts['h'][i]
            a = np.radians(self.box_extracts['a'][i])

            # get the box
            boxes.append(Polygon(get_box_edges(x, y, w, h, a)[:4]))

        return boxes

    def plot(self, ax, plot_sigma=True):
        '''
            Plot the the data for this jet object. Plots the
            start and end clustered points, and the associated
            extracts. Also plots the clustered and extracted
            box. Also plots a vector from the base to the top of the box

            Input
            -----
            ax : `matplotlib.Axes()`
                axis object to plot onto. The general use case for this
                function is to plot onto an existing axis which already has
                the subjet image and potentially other jet plots

            Outputs
            -------
            ims : list
                list of `matplotlib.Artist` objects that was created for this plot
        '''
        boxplot, = ax.plot(*self.box.exterior.xy, '-', color='white',
                             linewidth=0.8, zorder=10)
        startplot, = ax.plot(*self.start, 'bx', markersize=2, zorder=10)
        endplot, = ax.plot(*self.end, 'rx', markersize=2, zorder=10)

        start_ext = self.get_extract_starts()
        end_ext = self.get_extract_ends()

        # plot the extracts (start, end, box)
        startextplot, = ax.plot(
            start_ext[:, 0], start_ext[:, 1], 'k.', markersize=1.)
        endextplot, = ax.plot(
            end_ext[:, 0], end_ext[:, 1], 'k.', markersize=1.)
        boxextplots = []
        for box in self.get_extract_boxes():
            iou = box.intersection(self.box).area/box.union(self.box).area
            boxextplots.append(
                ax.plot(*box.exterior.xy, '-', color='limegreen', linewidth=0.5, alpha=0.65*iou+0.05)[0])

        # find the center of the box, so we can draw a vector through it
        center = np.mean(np.asarray(self.box.exterior.xy)[:, :4], axis=1)

        # create the rotation matrix to rotate a vector from solar north tos
        # the direction of the jet
        rotation = np.asarray([[np.cos(self.angle), -np.sin(self.angle)],
                               [np.sin(self.angle), np.cos(self.angle)]])

        # create a vector by choosing the top of the jet and base of the jet
        # as the two points
        point0 = center + np.matmul(rotation, np.asarray([0, self.height/2.]))
        point1 = center + np.matmul(rotation, np.asarray([0, -self.height/2.]))
        vec = point1 - point0

        base_points, height_points = self.get_width_height_pairs()

        arrowplot = ax.arrow(
            *point0, vec[0], vec[1], color='white', width=2, length_includes_head=True, head_width=10)

        if plot_sigma:
            if hasattr(self, 'sigma'):
                # calculate the bounding box for the cluster confidence
                plus_sigma, minus_sigma = sigma_shape(
                    self.cluster_values, self.sigma)

                # get the boxes edges
                plus_sigma_box = get_box_edges(*plus_sigma)
                minus_sigma_box = get_box_edges(*minus_sigma)

                # create a fill between the - and + sigma boxes
                x_p = plus_sigma_box[:, 0]
                y_p = plus_sigma_box[:, 1]
                x_m = minus_sigma_box[:, 0]
                y_m = minus_sigma_box[:, 1]
                ax.fill(
                    np.append(x_p, x_m[::-1]), np.append(y_p, y_m[::-1]), color='white', alpha=0.3)

        return [boxplot, startplot, endplot, startextplot, endextplot, *boxextplots, arrowplot]

    def autorotate(self):
        '''
            Find the rotation of the jet wrt to solar north and
            find the base width and height of the box
        '''
        box_points = np.transpose(self.box.exterior.xy)[:4, :]

        # find the distance between each point and the starting base
        dists = [np.linalg.norm((point - self.start)) for point in box_points]
        sorted_dists = np.argsort(dists)

        # the base points are the two points closest to the start
        base_points = np.array(
            [box_points[sorted_dists[0]], box_points[sorted_dists[1]]])

        # the height points are the next two
        rolled_points = np.delete(np.roll(box_points, -sorted_dists[0],
                                          axis=0), 0, axis=0)

        # we want to make sure that the order of the points
        # is in such a way that the point closest to the base
        # comes first -- this will ensure that the next point is
        # along the height line
        if np.linalg.norm(rolled_points[0, :]-base_points[1, :]) == 0:
            height_points = rolled_points[:2]
        else:
            height_points = rolled_points[::-1][:2]

        self.base_points = base_points
        self.height_points = height_points

        # also figure out the angle and additional size metadata
        # the angle is the angle between the height points and the base
        dh = self.height_points[1] - self.height_points[0]
        self.angle = np.arctan2(dh[0], -dh[1])

        self.height = np.linalg.norm(dh)
        self.width = np.linalg.norm(self.base_points[1] - self.base_points[0])

    def get_width_height_pairs(self):
        '''
            Outputs the base points and the height line segment
            points

            Outputs
            -------
            base_points : `numpy.ndarray`
                the pair of points that correspond to the base of the jet
            height_points : `numpy.ndarray`
                the pair of points that correspond to the height of the jet
        '''

        return self.base_points, self.height_points

class JetCluster:
    def __init__(self, jets):
        '''
            Initiate the JetCluster with a list of jet objects that are contained by that cluster.
        '''
        self.jets = jets

    def adding_new_attr(self, name_attr, value_attr):
        '''
            Add new attributes to the JetCluster
        Inputs
        ------
            name_attr: str
                name of the to be added property
            value_attr: any
                value of the to be added property
        '''
        setattr(self, name_attr, value_attr)

    def create_gif(self, output):
        '''
            Create a gif of the jet objects showing the
            image and the plots from the `Jet.plot()` method
        Inputs
        ------
            output: str
                name of the exported gif
        '''
        fig, ax = plt.subplots(1, 1, dpi=250)

        # create a temp plot so that we can get a size estimate
        subject0 = self.jets[0].subject

        ax.imshow(get_subject_image(subject0, 0))
        ax.axis('off')
        fig.tight_layout(pad=0)

        # loop through the frames and plot
        ims = []
        for jet in tqdm.tqdm(self.jets):
            subject = jet.subject
            for i in range(15):
                img = get_subject_image(subject, i)

                # first, plot the image
                im1 = ax.imshow(img)

                # for each jet, plot all the details
                # and add each plot artist to the list
                jetims = jet.plot(ax, plot_sigma=False)

                # combine all the plot artists together
                ims.append([im1, *jetims])

        # save the animation as a gif
        ani = animation.ArtistAnimation(fig, ims)
        ani.save(output, writer='imagemagick')

    def json_export(self, output):
        '''
            export one single jet cluster to output.json file
            Inputs
            ------
            output : str
                name of the exported json file
        '''
        json_export_list([self], output)

