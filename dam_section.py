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
        return (self._x, self._y)
    
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
    def get_points(self):
        """Returns the end points."""
        return self._points
        
    def get_end1(self):
        """Returns the first endpoint."""
        return self._points[0]
    
    def get_end2(self):
        """Returns the second endpoint."""
        return self._points[1]
        
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
        # Point of interest.
        xo, yo = point.get_xy()
        
        # Slope and intercept of line.
        m1 = self.get_slope()
        b1 = self.get_intercept()
        
        # Slope and intercept of orthogonal line to point.
        m2 = -1/m1
        b2 = yo - m2 * xo
        
        # Intersection of the two lines (projection point)
        x_intersec = (b2 - b1) / (m1 - m2)
        y_intersec = m1 * x_intersec + b1
        
        return  Line(point, Point(x_intersec, y_intersec))     
    
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
    """
    def __init__(self, name, data, elev, lat, long):
        """Initializes a Borehole instance."""
        self._name = name
        self._data = data
        self._elev = elev
        self._lat = lat
        self._long = long
    
    # Public getters.
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
    
    def get_point(self):
        """Returns the borehole coordinates as a point."""
        return Point(self._lat, self._long)
    
    # Public methods.
    def get_projection_point(self, line):
        """Returns the projection point of the borehole to a line of interest
        (i.e., orthogonal projection to a line of interest).
        
        Parameters:
            line : Line object representing the line of interest.
        """
        projection_line = line.get_orthogonal_line(self.get_point())
        return projection_line.get_intersection(line)
    
    def get_projection_line(self, line):
        """Returns a new Line object representing the orthogonal projection
        from the borehole location to a line of interest. 
        
        Parameters:
            line : Line object representing the line of interest.
        """
        return line.get_orthogonal_line(self.get_point())
    
    def get_section_loc(self, section_line, ref_point):
        """Returns the location of the borehole on a section.
        
        Parameters:
            section_line : Line object, representing the section line.
            ref_point : a reference point on the section_line from which to 
                to calculate the horizontal offset, distance or location from
                the reference. Should be an upstream or most distal point since
                we are deriving distances.
        """
        # The boring location.
        boring_point = self.get_point()
        
        # Get the location of the boring orthogonal projection to the section.
        projection_point = self.get_projection_point(section_line)
    
        # As ref_point and projection_point are both on the section_line, get
        # the horizontal distance between the two (in decimal degrees dd).
        x_dd = ref_point.get_distance(projection_point)
        
        # Convert from decimal degrees to ft. One degree of lattitude approx. 
        # equals 364,567.2 ft (69.1 miles).
        x_ft = x_dd * 364567.2 
        return x_ft
    
    def plot_borehole(self):
        """Plots the borehole."""
        pass


class BoringCollection:
    """BoringCollection representing a collection of individual boreholes of
    the same type.
    
    The idea is to plot collections of similar boreholes, e.g.
    mud rotary lithologies and classifications, or percent gravel/sand/fines,
    or continuous CPT data.
    
    Attributes:
        _borings : dictionary representing individual boreholes, where the keys
            are the unique identifiers.
        _boring_data : dataframe of raw boring data including lat/long.
        _boring_color_data : dataframe representing custom color mapping.
        _boring_interval_data : dataframe of raw boring interval data. 
        
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
    
    def get_color_data(self):
        """Returns the custom colormpa if it exists, otherwise sets and returns
        a generic colormap."""
        # If an existing colormap exists, return it.
        if self._boring_color_data:
            return self._boring_color_data
        # Otherwise, return a generic colormap.
        else:
            keys = self._interval_data.Classification.unique()
            cmap = cm.get_cmap(name="Spectral")
            cmap_segmented = cmap(np.linspace(0, 1, len(keys))) 
            colors = dict()
            for key, c in zip(keys, cmap_segmented):
                colors[key] = c
            self._boring_color_data = colors
            return self._boring_color_data

    # Public setters.
    def add_boring(self, boring_id, borehole):
        """Adds a boring to the collection."""
        self._borings[boring_id] = borehole
        
    def remove_boring(self, boring_id):
        """Removes a boring from the collection if it exists."""
        if self._borings.get(boring_id, None):
            del self._borings[boring_id]

    def set_color_data(self, df_color):
        """Sets custom color data.
        Parameters:
            df_color : dataframe of classification/lithology to color.
        """
        # Convert the color data from dataframe to dict with classification
        # values as keys.
        color_data = df_color.set_index('Classification')['Color'].to_dict()
        self._boring_color_data = color_data
    
    def remove_color_data(self):
        """Resets/removes custom color data."""
        self._boring_color_data = None
    
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
        
        self._boring_data = boring_data
        self._interval_data = interval_data


class Section:
    """Section class representing dam section.
    
    Can be used to plot a dam section, and overlay/project boring data from
    BoringCollections.
    
    Attributes:
        _data : nested dictionary representing the section data where the 
                top-level keys are 'sections' and 'colors', referring to the
                geometry and color data, respectively. The bottom-level keys
                correspond to the geometry section names/columns.
                referring 
        _line : Line object representing the projection axis/line/dam 
                alignment.
        _collections : dictionary representing various BoreholeCollection
                objects, where the key corresponds to a unique identifer.
    """
    def __init__(self):
        self._data = {'sections': {}, 'colors':{}}
        self._line = None
        self._collections = {}
    
    # Public getters.
    def get_data(self):
        """"Returns the section and color data dict."""
        return self._data
    
    def get_line(self):
        """Returns the line object of the Section."""
        return self._line
    
    def get_zones(self):
        """Returns the keys/zones of the section data."""
        return self._data['sections'].keys()
    
    def get_section(self, key):
        """Returns the section data for the specified key."""
        return self._data['sections'].get(key, None)
    
    def get_color(self, key):
        """Returns the section color for the specified key"""
        return self._data['colors'].get(key, None)
    
    # Public setters. 
    def set_line(self, line):
        """Sets the alignment of the section from a Line object."""
        self._line = line
        
    def set_line_from_points(self, p1, p2):
        """Sets the alignment of the section from two Point objects."""
        self._line = Line(p1, p2)
        
    def add_collection(self, collection_id, boring_collection):
        """Adds a boring collection.
        
        Parameters:
            collection_id : unique identifier of the borehole collection. 
            boring_collection : BoringCollection object representing the group
                of borings to add. 
        """
        self._collections[collection_id] = boring_collection
        
    def remove_collection(self, collection_id):
        """Removes a boring collection.
        
        Parameters:
            collection_id : unique identifier of the borehole collection. 
        """
        if self._collections.get(collection_id, None):
            del self._collections[collection_id]
            
    def return_boring_collection(self, collection_id):
        """Returns a boring collection object given the unique identifier.
        
        Parameters:
            collection_id : unique identifier of the borehole collection.
        """
        return self._collections.get(collection_id, None)
    
    # Public methods.
    def read_section_data(self, fname, sheet_name=None):
        """Reads dam section data from a csv or excel file.
        
        Parameters:
            fname : string representing the filename.
            sheet_name : string representing sheet name, if excel file.
        """
        try: 
            ext = fname.split('.')[-1]
            if ext in ('.xlsx','.xlsm', '.xlsb', '.xltx'):
                table = pd.read_excel(fname, sheet_name)
            else:
                table = pd.read_csv(fname)
            
            columns = table.columns
            # Reads every 2 columns as a set until we reach the end of the table.
            for i in range(0, len(columns), 2):
                
                # Gets the column name/key.
                key = columns[i]
                
                # Extracts the geometry tag, and color assginment.
                tag = table.iloc[0,i+1]
                color = table.iloc[1,i+1]
                
                # Gets the table x and y labels.
                x_lab = table.iloc[2,i]
                y_lab = table.iloc[2,i+1]
                
                # Gets the x and y data.
                x = table.iloc[3:,i]
                y = table.iloc[3:,i+1]
                
                # Creates a new dataframe of x, y data, striping all NAN rows, 
                # and casting to float type.
                df = pd.DataFrame(data={x_lab: x, y_lab: y})
                df = df.dropna(axis=0, how='all').astype(float)
        
                # Assigns the section dataframes and colors by section_name to 
                # attributes.
                self._data['sections'][key] = df
                self._data['colors'][key] = color
        except:
            print("Failed to load section data.")
                
    
    def plot_sections(self, ax=None, fig=None, x_scale=1, y_scale=1):
        """Plots the section data, and returns the figure and axis handle.
        
        Certain keywords (case-insensitive) are resevered, e.g.
        'profile' or 'dam profile' - for the dam profile is plotted as a 
            black line.
        'resevoir level', or 'water level' - for the water level is plotted as
            a blue line.
            
        Parameters:
            ax  : matplotlib axis handle.
            fig : matplotlib figure handle.
            x_scale : scales the horizontal geometry by x_scale.
            y_scale : scales the vertical geometry by y_scale.
        """
        try:
            # Sets up a new figure and axis handle if we aren't given one.
            if ax == None:
                fig, ax = plt.subplots()
                ax.set_aspect('equal','box')
            
            # Plots all of the section data geometries.
            for key in self.get_zones():
                
                section = self.get_section(key)
                color = self.get_color(key)
                
                if key.lower() in ('dam profile', 'profile'):
                    ax.plot(section.x * x_scale, section.elev * y_scale,
                            c='k', lw=1.0, zorder=0, label=key)
                elif key.lower() in ('resevoir level', 'water level'):
                    ax.plot(section.x * x_scale, section.elev * y_scale,
                            c='b', lw=1.0, zorder=0, label=key)
                else:
                    ax.fill(section.x * x_scale, section.elev * y_scale,
                            c=color,  lw=1.0, zorder=0, label=key)
            return fig, ax
        except:
            print("Error when plotting dam section.")
            
            
    def add_section_legend(self, ax, 
                           loc='upper center', 
                           fontsize=10, 
                           strip=()):
        """Adds a legend for sections to an axes.
        
        Parameters:
            ax : matplotlib axes to add legend to.
            loc : legend location.
            fontsize : legend font size.
            strip : tuple/list of strings segments to strip/delete from legend
                labels.
        """
        # Get the handles and legend labels
        handles, labels = ax.get_legend_handles_labels()
        
        # Delete/remove specific string segments
        modified_labels = []
        for lab in labels:
            for string in strip:
                lab = lab.replace(string,'')
            lab = lab.rstrip()
            modified_labels.append(lab)
        
        # Adds legend with repeating labels 
        by_label = dict(zip(modified_labels, handles))
        ax.legend(by_label.values(), by_label.keys(),
                  ncol=len(by_label), loc=loc, fontsize=fontsize)
        return ax
                 
    
    def plot_borings(self, collection_id, ax=None, fig=None, 
                     boring_offset = 0,
                     boring_width = 10,
                     bh_font_size = 5,
                     bh_label_voffset = 5,
                     ):
        """Plots the boring data.
        
        Parameters:
            collection_id : unique identifier representing boring collection.
            ax : passed matplotlib axes handle.
            fig : passed matplotlib figure handle.
            boring_offset : specifies a custom horizontal offset/shift applied
                to all boreholes, in ft.
            boring_width : width of the boreholes, in ft.
            bh_font_size : font size of borehole descriptions.
            bh_label_voffset : vertical shift of borehole id, in ft.
        """
        # Get the boring data and colors
        boring_collection = self._collections[collection_id]
        borings = boring_collection.get_borings()
        colors = boring_collection.get_color_data()
        
        try:
            # Sets up a new figure and axis handle if we aren't given one.
            if ax == None:
                fig, ax = plt.subplots()
                ax.set_aspect('equal','box')
            
            # Loop over the borings.
            for bh in borings:
                
                # Get the boring data and location/horizontal offset.
                df = borings[bh].get_data()
                xo = borings[bh].get_section_loc(self._line, 
                                                 self._line.get_end1())
                
                # Apply custom offset to borings. 
                xo += boring_offset
                
                # Add boring label.
                ax.text(xo, df.iloc[0]['Start Elev.'] + bh_label_voffset, bh)
                
                # Plot the interval data.
                for i, row in df.iterrows():
                    
                    x = np.array([xo, xo + boring_width]).flatten()
                    y1 = [row['End Elev.'],row['End Elev.']]
                    y2 = [row['Start Elev.'],row['Start Elev.']]
                    col = colors.get(row.Classification, 'w')
                    ax.fill_between(x, y1, y2=y2, color=col, ec='k', lw=0.5)
                    
                    # Add a borehole description
                    try:        
                        description = row.Description.split(",")[0]
                    except:
                        description = ""
                    ax.text(xo + boring_width,
                            (row['End Elev.'] + row['Start Elev.'])/2, 
                             f"{row.Classification} ({description})", 
                             fontsize = bh_font_size
                             )
            return fig, ax
        except:
            print("Error when plotting borehole collection.")
            

def main():
    # Boring inputs.
    boring_ids = [5, 10, 23, 29, 22, 11, 1, 3, 2, 4]
    fname = "Terminal Dam Boring Intervals.xlsx"
    boring_data = pd.read_excel(fname, sheet_name="borings")
    interval_data = pd.read_excel(fname, sheet_name="boring-data")
    boring_color_data = pd.read_excel(fname, sheet_name="colors")
    
    # Create the boring collection and add borings, and color data.
    bhc = BoringCollection()
    bhc.process_borings(boring_ids, boring_data, interval_data)
    bhc.set_color_data(boring_color_data)
    
    # Create the dam section, and define alignment.
    section = Section()
    section.read_section_data("Terminal Dam Section.csv")
    us_toe = Point(35.1708928, -120.5345804)
    ds_toe = Point(35.1696816, -120.5333135)
    section.set_line_from_points(us_toe, ds_toe)
    
    # Attach the boring collection to the dam section.
    section.add_collection('wcs-borings', bhc)
    
    # Plot the dam section and borings.
    fig, ax = plt.subplots(figsize=(15,15))
    section.plot_sections(ax=ax)
    section.plot_borings('wcs-borings', ax=ax, boring_offset=-15)
    ax.set_ylim(bottom=150, top=380)
    ax.set_xlim(left=-25, right=850)
    ax.set_aspect('equal','box')    
    ax.set_xlabel('Distance (ft)')
    ax.set_ylabel('Elevation (ft)')
    section.add_section_legend(ax, strip=('US', 'DS'))
    fig.savefig("test.svg", dpi=300, bbox_inches="tight", format="svg")
    plt.show()
    

if __name__ == '__main__':
    main()
    