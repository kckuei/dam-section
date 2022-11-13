"""
By: kckuei

Program Description:  Implements several classes for plotting/visualizing 
borehole data on dam sections. The following classes are implemented including
Point, Line, Borehole, BoringCollection, and Section classes.

The idea is that discrete continuous or interval data can be represented by
Borehole objects. Borehole object locations can be represented by Point
objects. Dam cross-sections can be represented by a Section object. A section
has an inherent direction or alignment which can be represented by a Line
object, and in turn Point objects.

We can group collections of Borehole objects as a BoringCollection, and
assign them to Section objects for plotting. To plot a Borehole object or
BoringCollection object, we project them (from their point) to the section
(the Line).

"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm



class Point:
    """Point class for representing coordinates of boreholes. Points can 
    also be used to define the line or alignment that a Section is define don.
    
    Attributes:
        _x : x coordinate value (in plan view lattitude).
        _y : y-coordinate value (in plan view longitude).
    """
    def __init__(self, x, y):
        """Initializes a Point object."""
        self._x = x
        self._y = y
    
    # Public getters.
    def get_x(self):
        """Returns x value."""
        return self._x
    
    def get_y(self):
        """Returns y value."""
        return self._y
    
    def get_xy(self):
        """Returns (x,y) tuple."""
        return (self._x, self.y_)
    
    # Public methods.
    def get_distance(self, p):
        """Get distance between current and provided Point objects using
        distance formula. 
        
        Parameters:
            p : Point object.
        """
        x1, y1 = self.get_xy()
        x2, y2 = p.get_xy()
        dist = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
        return dist
    
    def get_slope_intercept(self, p):
        """Returns the slope and intercept as a tuple of the current point
        with the supplied point.
        
        Parameters:
            p : point object.
        """
        x1, y1 = self.get_xy()
        x2, y2 = p.get_xy()
        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1
        return slope, intercept 
    

class Line:
    """Line class for representing the direction/alignment of a section. 
    
    Points can be projected to a line.
    
    Attributes:
        _points     : list of Point objects representing the seed points for
                        the line.
        _slope      : slope of the line (in plan view lat/long).
        _intercept  : intercept of the line (in plan view lat/long).
    """
    def __init__(self, p1, p2):
        """Initializes a Line object.""" 
        slope, intercept = p1.get_slope_intercept(p2)
        self._points = [p1, p2]
        self._slope = slope
        self._intercept = intercept
      
    # Puclic getters.
    def get_slope(self):
        """"Returns the slope."""
        return self._slope
      
    def get_intercept(self):
        """Returns the intercept."""
        return self._intercept
      
    # Public methods.
    def feval(self, x):
        """Evaluates f(x) from the slope and intercept of the line.
        
        Parameters:
            x : horizontal offset/distance along the line.
        """
        return self._slope * x + self._intercept

    def get_orthogonal_line(self, point):
        """Returns a new Line object representing the orthogonal projection 
        from the current line to the passed Point object.
        
        Parameters:
            point : Point object.
        """
        m1 = self.get_slope()
        m2 = -1/m1
        xo, yo = point.get_xy()
        b2 = yo - m2 * xo
        return Line(m2, b2)
    
    def get_intersection(self, line):
        """Returns the intersection of the current line object with the 
        passed Line object as a Point object.
        
        Parameters:
            line : Line object.
        """
        m1, b1 = self.get_slope(), self.get_intercept()
        m2, b2 = line.get_slope(), line.get_intercept()
        x_intersec = (b2 - b1) / (m1 - m2)
        y_intersec = m1 * x_intersec + b1
        return Point(x_intersec, y_intersec)
      
    # Static method
    def get_line_from_points(p1, p2):
        """Returns a new Line object, given two Point objects.
        
        Parameters:
            p1 : Point object 1. 
            p2 : Point object 2.
        """
        return Line(*p1.get_slope_intercept(p2))



class Borehole:
    """Borehole class for representing different kinds of soundings, e.g. CPT, 
    borehole, etc. 
    
    Idea is to allow plotting continuous/discrete or interval-like data such as
    borehole lithology eventually.
    
    Attributes:
        _name    : unique name or id
        _data    : borehole data
        _elev    : borehole elevation datum
        _lat     : borehole latittude 
        _long    : borehole longitude
        
        _projection_point : projection point to a particular section.          ; Maybe should not be attributes of the borehole.
        _projection_line  : projection line to a particular section.           ; Could be methods of a section given a boring. 
        _section_loc      : projection x-loc on the particular section.        ; Borehole can be attached to multiple sections. 
    """
    def __init__(self, name, data, elev, lat, long):
        """Initializes a Borehole instance."""
        self._name = name
        self._data = data
        self._elev = elev
        self._lat = lat
        self._long = long
        self._projection_point = None
        self._projection_line = None
        self._section_loc = None
    
    # Public getters
    def get_data(self): 
        """"Returns the data."""
        return self._data
    
    def get_elev(self):
        """Returns the elevation."""
        return self._elev
    
    def get_lat(self):
        """Returns the lattiude."""
        return self._lat
    
    def get_long(self):
        """Returns the longitude."""
        return self._long
    
    def get_projection_point(self):
        """Returns the projection point."""
        return self._projection_point
    
    def get_projection_line(self):
        """Returns the projection line."""
        return self._projection_line
    
    def get_section_loc(self):
        """Returns the location on the section."""
        return self._section_loc
    
    # Public setters
    def set_projection_point(self, point):
        """Sets the projection point to a section."""
        self._projection_point = point
        
    def set_projection_line(self, line):
        """Sets the projection line."""
        self._projection_line = line
    
    def set_section_loc(self, x):
        """Sets the location on the section."""
        self._section_loc = x



class BoringCollection:
    """BoringCollection representing a collection of individual boreholes of
    the same type.
    
    The idea is to plot collections of similar boreholes, e.g.
    mud rotary lithologies and classifications, or percent gravel/sand/fines,
    or continuous CPT data.
    
    Attributes:
        _borings :
        _boring_data :
        _boring_color_data :
        _boring_interval_data :
        
    """
    def __init__(self):
        self._borings = {}
        self._boring_data = None
        self._boring_color_data = None
        self._interval_data = None
        
    # Public getters.
    def get_borings(self):
        """Returns the borings."""
        return self._borings
    
    def get_boring_ids(self):
        """Returns the boring ids."""
        return self._borings.keys()
    
    # Public methods.
    def process_borings(self, boring_ids, boring_data, interval_data):
        """Process the borings, and return a dictionary of borings.
        Parameters:
            boring_ids : list of integers representing borings to process.
            boring_data : Dataframe representing boring data. 
            interval_data : Dataframe representing interval data. 
        """
        for boring_id in boring_ids:
            
            data = interval_data[interval_data.Hole == boring_id].copy(deep=True)
            
            hole_meta = boring_data[boring_data.Hole == boring_id]
            
            elev = hole_meta['Ground Elev.'].values
            lat = hole_meta['Lattitude'].values
            long = hole_meta['Longitude'].values
            
            data['Start Elev.'] = elev - data['Start Depth']
            data['End Elev.'] = elev - data['End Depth']
            
            self._borings[boring_id] = Borehole(boring_id, data, elev, lat, long)

    
    def update_boring_projections(borings, section_line, xy_us_toe):
        """Updates boring projections."""  
        # Add the boring projections to the boring objects
        for bh in borings:
            
            # Get the boring location
            xy = (borings[bh].get_lat(), borings[bh].get_long())
        
            # Get borehole projection line by projecting from borehole location to section line
            projection_line = get_orthogonal_line(section_line, xy)
            
            # Get the intersection point, or borehole location projected to section
            projected_point = get_intersection(section_line, projection_line)
            
            # Get the in-plane distance to the borehole on the section
            x = get_distance(xy_us_toe, projected_point)
        
            # One degree of lattitude equals ~ 364,567.2 ft (69.1 miles); approximation; revisit later and check
            x_ft = x * 364567.2 
        
            borings[bh].set_section_loc(x_ft)
            borings[bh].set_projection_point(projected_point)
            borings[bh].set_projection_line(projection_line)


class Section:
    """Section class representing any dam section (e.g. transverse,
    longitudinal, etc.).
    """
    def __init__(self, line):
        self._line = line
        self._section_data = None
        self._color_data = None
        self._boringcollection = None
    
    def add_section_data(self, excel_file, sheet_name):
        """Adds dam section datafrom an excel file and returns a dictionary 
        with each section represented by a dataframe.
        """
        
        raw = pd.read_excel(excel_file, sheet_name=sheet_name)
        
        section_data = {}
        color_data = {}
        
        columns = raw.columns
        for i in range(0, len(columns), 2):
            section_name = columns[i]
            
            tag = raw.iloc[0,i+1]
            color = raw.iloc[1,i+1]
            
            x_lab = raw.iloc[2,i]
            y_lab = raw.iloc[2,i+1]
            
            x = raw.iloc[3:,i]
            y = raw.iloc[3:,i+1]
            
            df = pd.DataFrame(data={x_lab: x, y_lab: y}).dropna()
    
            section_data[section_name] = df
            color_data[section_name] = color
            
        self._section_data = section_data
        self._color_data = color_data
        
        return self._section_data, self._color_data

    def plot_sections(fig, ax, section_data, color_data):
        """Plots the section data."""
        
        for key in section_data.keys():
            
            section = section_data[key]
            
            if 'profile' in key.lower():
                ax.plot(section.x, section.elev,
                        c='k', lw=1.0, zorder=0, label=key)
            elif 'resevoir level' in key.lower():
                ax.plot(section.x, section.elev,
                        c='b', lw=1.0, zorder=0, label=key)
            else:
                ax.fill(section.x, section.elev,
                        c=color_data[key], lw=1.0, zorder=0, label=key)
                
        return fig, ax
    
    
    def plot_borings(fig, ax, borings, colors):
        """Plots the boring data."""
        width = 10
        fs = 5
        bh_label_voffset = 5
        
        for bh in borings:
            df = borings[bh].get_data()
            xo = borings[bh].get_section_loc()
            plt.text(xo, df.iloc[0]['Start Elev.'] + bh_label_voffset, bh)
            for i, row in df.iterrows():
                x = np.array([xo, xo + width]).flatten()
                y1 = [row['End Elev.'],row['End Elev.']]
                y2 = [row['Start Elev.'],row['Start Elev.']]
                plt.fill_between(x, y1, y2=y2, color=colors[row.Classification],
                                 ec='k', lw=0.5)
                try:        
                    description = row.Description.split(",")[0]
                except:
                    description = ""
                plt.text(xo + width, (row['End Elev.'] + row['Start Elev.'])/2, 
                         f"{row.Classification} ({description})", 
                         fontsize = fs)


    def check_projection(borings, xy_us_toe, xy_ds_toe):
        """Checks the boring projection."""
        fig, ax = plt.subplots()
        for bh in borings:
            plt.plot(borings[bh].get_lat(), borings[bh].get_long(), 'o')
        for bh in borings:
            x, y = borings[bh].get_projection_point()
            plt.plot(x, y, 'o', mfc='w')
        xlim = plt.gca().get_xlim()
        for bh in borings:
            plt.plot(xlim,
                [borings[bh].get_projection_line().f(xlim[0]), 
                 borings[bh].get_projection_line().f(xlim[1])]) 
        x, y = list(zip(*[xy_us_toe, xy_ds_toe]))
        plt.plot(x,y)
        plt.axis('equal')

    
    def format_axes(fig, ax):
        """"Format axes."""
        ax.set_aspect('equal','box')
        # ax.set_aspect(2.0)
        fig.tight_layout()
        # ax.set_ylim(bottom=100, top=360)
        ax.set_ylim(bottom=150, top=360)
        ax.set_xlim(left=-25, right=800)
        ax.set_xlabel('Distance (ft)')
        ax.set_ylabel('Elevation (ft)')
        
        # Legend with repeating labels removed (ignoring US and DS tags)
        handles, labels = ax.get_legend_handles_labels()
        labels = [lab.replace('US','').replace('DS','').rstrip() for lab in labels]
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), ncol=len(by_label), loc='upper center')



def main():
    
    # Load the data
    fname = "Terminal Dam Boring Intervals.xlsx"
    boring_data = pd.read_excel(fname, sheet_name="borings")
    interval_data = pd.read_excel(fname, sheet_name="boring-data")
    boring_color_data = pd.read_excel(fname, sheet_name="colors")
    boring_ids = [5, 10, 23, 29, 22, 11, 1, 3, 2, 4]
    
    
    bhc = BoringCollection()
    bhc.process_borings(boring_ids, boring_data, interval_data)
    
    
    print('hello')
    
    # # Holes to process
    # boring_ids = [5, 10, 23, 29, 22, 11, 1, 3, 2, 4]
    
    # # Process the borings
    # borings = process_borings(boring_ids, boring_data, interval_data)
    
    # # Get section line
    # # xy_us_toe = (35.1708928, -120.5345804)
    # # xy_ds_toe = (35.1696816, -120.5333135)
    # xy_us_toe = (35.1707027, -120.5346784)
    # xy_ds_toe = (35.1695906, -120.5335153)
    # section_line = get_line_from_points(xy_us_toe, xy_ds_toe)
        
    # # Update boring projections
    # update_boring_projections(borings, section_line, xy_us_toe)
    
    # # Create colormap
    # ## sorted(list(filter(lambda x: type(x) is str, key_values)))
    # # key_values = interval_data.Classification.unique()
    # # cmap = cm.get_cmap(name="Spectral")
    # # cmap_segmented = cmap(np.linspace(0, 1, len(key_values))) 
    # # colors = dict()
    # # for key, c in zip(key_values, cmap_segmented):
    # #     colors[key] = c
    # colors = boring_color_data.set_index('Classification')['Color'].to_dict()
    
    # # Setup figure
    # fig, ax = plt.subplots(figsize=(15,15))
    
    # # Get the section data
    # section_data, color_data = read_section_data("Terminal Dam Section.xlsx", "Long")
    
    # # Draw the dam section, etc.
    # plot_sections(fig, ax, section_data, color_data)
    # plot_borings(fig, ax, borings, colors) 
    # format_axes(fig, ax)
    
    # # Save figure
    # fig.savefig("test.svg", dpi=300, bbox_inches="tight", format="svg")
    # plt.show()


if __name__ == '__main__':
    main()
    