import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm


class Borehole:
    """Basic borehole class."""
    def __init__(self, name, data, elev, lat, long):
        self._name = name
        self._data = data
        self._elev = elev
        self._lat = lat
        self._long = long
        self._projection_point = None
        self._projection_line = None
        self._section_loc = None
    
    def get_data(self): 
        return self._data
    
    def get_elev(self):
        return self._elev
    
    def get_lat(self):
        return self._lat
    
    def get_long(self):
        return self._long
    
    def get_projection_point(self):
        return self._projection_point
    
    def get_projection_line(self):
        return self._projection_line
    
    def get_section_loc(self):
        return self._section_loc
    
    def set_projection_point(self, point):
        self._projection_point = point
        
    def set_projection_line(self, line):
        self._projection_line = line
    
    def set_section_loc(self, x):
        self._section_loc = x


class Line:
  """Basic line class."""
  def __init__(self, slope, intercept):
    self._slope = slope
    self._intercept = intercept
    
  def f(self, x):
    return self._slope * x + self._intercept
    
  def get_slope(self):
    return self._slope
    
  def get_intercept(self):
    return self._intercept


def get_line_from_points(p1, p2):
  """Get line from two points p1 and p2."""
  x1, y1 = p1
  x2, y2 = p2
  m = (y2 - y1) / (x2 - x1)
  b = y1 - m * x1
  return Line(m, b)


def get_orthogonal_line(line, point):
  """Gets line orthogonal to slope_intercept and passing through point."""
  m1 = line.get_slope()
  m2 = -1/m1
  xo, yo = point
  b2 = yo - m2 * xo
  return Line(m2, b2)


def get_intersection(line1, line2):
  """Get intersection point of two lines."""
  m1, b1 = line1.get_slope(), line1.get_intercept()
  m2, b2 = line2.get_slope(), line2.get_intercept()
  x_intersec = (b2 - b1) / (m1 - m2)
  y_intersec = m1 * x_intersec + b1
  return (x_intersec, y_intersec)
  

def get_distance(p1, p2):
  """Get distance between two points p1 and p2."""
  x1, y1 = p1
  x2, y2 = p2
  dist = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
  return dist





def read_section_data(excel_file, sheet_name):
    """Reads dam section data and returns a dictionary with 
    each section represented by a dataframe."""
    
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
    return section_data, color_data


def process_borings(boring_ids, boring_data, interval_data):
    """Process the borings, and return a dictionary of borings."""
    borings = {}
    for boring_id in boring_ids:
        data = interval_data[interval_data.Hole == boring_id].copy(deep=True)
        hole_meta = boring_data[boring_data.Hole == boring_id]
        elev = hole_meta['Ground Elev.'].values
        lat = hole_meta['Lattitude'].values
        long = hole_meta['Longitude'].values
        data['Start Elev.'] = elev - data['Start Depth']
        data['End Elev.'] = elev - data['End Depth']
        borings[boring_id] = Borehole(boring_id, data, elev, lat, long)
    return borings 


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





def main():
    
    # Load the data
    fname = "Terminal Dam Boring Intervals.xlsx"
    boring_data = pd.read_excel(fname, sheet_name="borings")
    interval_data = pd.read_excel(fname, sheet_name="boring-data")
    boring_color_data = pd.read_excel(fname, sheet_name="colors")
    
    # Holes to process
    boring_ids = [5, 10, 23, 29, 22, 11, 1, 3, 2, 4]
    
    # Process the borings
    borings = process_borings(boring_ids, boring_data, interval_data)
    
    # Get section line
    # xy_us_toe = (35.1708928, -120.5345804)
    # xy_ds_toe = (35.1696816, -120.5333135)
    xy_us_toe = (35.1707027, -120.5346784)
    xy_ds_toe = (35.1695906, -120.5335153)
    section_line = get_line_from_points(xy_us_toe, xy_ds_toe)
        
    # Update boring projections
    update_boring_projections(borings, section_line, xy_us_toe)
    
    # Create colormap
    ## sorted(list(filter(lambda x: type(x) is str, key_values)))
    # key_values = interval_data.Classification.unique()
    # cmap = cm.get_cmap(name="Spectral")
    # cmap_segmented = cmap(np.linspace(0, 1, len(key_values))) 
    # colors = dict()
    # for key, c in zip(key_values, cmap_segmented):
    #     colors[key] = c
    colors = boring_color_data.set_index('Classification')['Color'].to_dict()
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(15,15))
    
    # Get the section data
    section_data, color_data = read_section_data("Terminal Dam Section.xlsx", "Long")
    
    # Draw the dam section, etc.
    plot_sections(fig, ax, section_data, color_data)
    plot_borings(fig, ax, borings, colors) 
    format_axes(fig, ax)
    
    # Save figure
    fig.savefig("test.svg", dpi=300, bbox_inches="tight", format="svg")
    plt.show()


if __name__ == '__main__':
    main()
    