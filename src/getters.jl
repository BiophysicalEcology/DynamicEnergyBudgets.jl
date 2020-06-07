# We use getters so the model never accesses struct fields directly, 
# except in methods for specific leaf types where declaring all the gettier 
# methods would be overkill.
#
# A benefit of this is parameters can be easily moved between shared and specific 
# parameter structs.
