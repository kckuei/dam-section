# Classes

# Points
#    The most rudimentary object  
# Lines
#    Lines are made from points
#    Borings can have a projection with respect to a line
# Borings
#    Can represent different kinds of soundings
#    A boring contains continuous or interval data
#    A boring has a location (point)
#    Has an elevation datum
# BoringCollection
#    A boring collection has a defined colormap
#    Has attributes for how it gets plotted
# Sections
#    Can be a plan, transverse, or longitudinal section
#    Defined by a collections of sections
#    The sections have a defined colormap
#    Can add or remove sections
#    A section has a defined line
#    A section can add or remove borings
#    Can add or remove sections
#    We can plot sections
#    Has attributes for how it gets plotted
# Figure/Canvas
#    We can add a dam section or a boring collection
#    Can have multiple boring collections, but only a single dam section to plot
#    We can format the axes
#    We can save it when we're done
#    Plots data differently if it is a plan vs transverse/longitudinal section